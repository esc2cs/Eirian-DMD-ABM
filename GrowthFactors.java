/**
 * 
 */
package muscle;

import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;

/**
 * @author Kelley Virgilio This class contains the growth factors for the model:
 *         getGF, setGF functions defined for all
 */
public class GrowthFactors {

	private static Grid<Object> grid;
	public static ContinuousSpace<Object> space;
	private static BufferedGridValueLayer mcpSpatial;
	public static int activeDelay = 96; // time delay from latent to active tgf, 96 normal
	public static double totalMCP;
	
	// DISEASE STATE PARAMETERS
	public static double m1MacAdded = 0; // 0 at healthy
	public static double m2MacAdded = 0; // 0 at healthy
	public static double mdxTGF = 0; // 0 at healthy
	// private static int activeDelay = 1; // void at healthy; add baseline active
	// TGF
	public static double inflamWeight = 0; // 0 at healthy 

	static final int numGrowthFactors = 32;// = 31; (w/o spatial MCP) // Number of growth factors in model
	public static double[] growthFactors = new double[numGrowthFactors]; // growth factors
	static final int simLength = 10000 + activeDelay + 1;
	public static double[] activeTgf = new double[simLength]; // holds active TGF at each time step
	
	public GrowthFactors(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, ContinuousSpace<Object> space) {
		this.grid = grid;
		this.space = space;
		this.mcpSpatial = mcpSpatial;
	}

	public static double[] getGrowthFactors() {
		return growthFactors;
	}

	// Growth factor solver
	public static double[] growthFactorSolver(double[] inflamCells, Context<Object> context) {
		// inflammCells: 0 RM; 1 N; 2 Na; 3 M1; 4 M1ae; 5 M1de; 6 M2
		double[] growthFactorsTemp = new double[numGrowthFactors];
		List<Object> activeFibroblasts = Fibroblast.getActiveFibroblasts(context);// active fibroblasts
		double numActiveFibrob = activeFibroblasts.size(); // number of active fibroblasts
		double inflamFibroblasts = 0; // fibroblasts only secrete certain factors in an inflammatory environment
		double inflamSSC = 0; // sscs only secrete certain factors in an inflammatory environment
		double numActiveSSC = SSC.getActiveSSCs(context).size(); // number of active sscs
		// exclude signals from fully differentiated SSCs
		double numActiveSecretingSSC = SSC.getNumActSecretingSSCs(context); // secreting sscs
		double numMyofbs = Fibroblast.getMyofibroblasts(context).size(); // number of myofibroblasts
		double resMac = 0.; // number of resident macrophages scaled to the amount of muscle damage
		double percentNecrotic = Necrosis.getPercentNecrotic(context); // percent of muscle that is necrotic
		// activation is a function of initial damage only

		//get total spatial MCP
		List<Object> ecms = ECM.getECM(context);
		List<Object> necrotic = Necrosis.getNecrosis(context);
		List<Object> fibers = Fiber.getFiberElems(context);
		
		double totalMCP = 0;
		for (Object obj : ecms) {
			GridPoint pt = grid.getLocation(obj);
			totalMCP = totalMCP + mcpSpatial.get(pt.getX(), pt.getY());
		}
		for (Object obj : necrotic) {
			GridPoint pt = grid.getLocation(obj);
			totalMCP = totalMCP + mcpSpatial.get(pt.getX(), pt.getY());
		}
		for (Object obj : fibers) {
			GridPoint pt = grid.getLocation(obj);
			totalMCP = totalMCP + mcpSpatial.get(pt.getX(), pt.getY());
		}
		growthFactors[31] = totalMCP;
		
		int timestep = 1;
		if (percentNecrotic < 0.001) { // lower limit of damage for resident macrophages to detect
			resMac = 0;
		} else {
			//resMac = inflamCells[0]; // else all resident macrophages detect damage
			resMac = Macrophage.getMres(context).size(); // get number of resMs
		}
		// INFLAMMATION WEIGHTING FUNCTION:
		// weighting function to determine if it is a pro-inflammatory or
		// anti-inflammatory environment
		inflamWeight = (growthFactors[1] - getActiveTgf(InflamCell.getTick())) / (40. * Fiber.origFiberNumber);
		if (inflamWeight < 0) {
			inflamWeight = 0;
		}
		inflamFibroblasts = numActiveFibrob * inflamWeight; // number of fibroblasts weighted by inflammatory
															// environment
		inflamSSC = numActiveSSC * inflamWeight; // number of sscs weighted by inflammatory environment

		//double m1Mac = inflamCells[3]; // m1 macrophages
		//double m2Mac = inflamCells[6]; // m2 macrophages
//		double m1Mac = Macrophage.getM1(context).size()/10; // m1 macrophages
//		double m2Mac = Macrophage.getM2(context).size()/10; // m2 macrophages
//		double m1ae = Macrophage.getM1ae(context).size()/10; // m1 macrophages that phagocytosed apoptotic neutrophils
//		double m1de = Macrophage.getM1ae(context).size()/10; // m1 macrophages that phagocytosed necrosis
//		double n = Neutrophil.getNeutrophils(context).size();
//		double na = Neutrophil.getApoptosed(context).size();


		// DISEASE STATE PARAMETERS
		// MACROPHAGE PHENOTYPE ANALYSIS:
		// M1: INFLAMMATORY-- add in a constant level of inflammatory M1s
		// m1Mac = m1Mac + m1MacAdded; // 0 added at healthy
		//
		// M2: anti-inflammatory- add in constant level of M2
		// m2Mac = m2Mac + m2MacAdded; // 0 added at healthy

		// LIT REPLICATION SCALING PARAMTERS: Murphy et al. 2011
		numActiveFibrob = numActiveFibrob * (1 / Fiber.tcf4Scale); // scale secretions --> less // tcf4scale = 1 at
																	// healthy = no change
		inflamFibroblasts = inflamFibroblasts * (1 / Fiber.tcf4Scale); // scale secretions --> less // tcf4scale = 1 at
																		// healthy = no change
		numActiveSecretingSSC = numActiveSecretingSSC * (1 / Fiber.pax7Scale); // scale secretions --> more // pax7scale
																				// = 1 at healthy = no change
		inflamSSC = inflamSSC * (1 / Fiber.pax7Scale); // scale secretions --> more // pax7scale = 1 at healthy = no
														// change
		numMyofbs = numMyofbs * (1 / Fiber.tcf4Scale); // scale secretions --> less // tcf4scale = 1 at healthy = no
														// change

		// INFLAMMATORY CELLS:
		// 0 RM- resident macrophages
		// 1 N- neutrophils
		// 2 Na- apoptotic neutrophils
		// 3 M1- M1 macrophages
		// 4 M1ae- M1 apoptotic eating
		// 5 M1de- M1 debris eating
		// 6 M2- M2 macrophages
	
		
		// Inflam cell agent counts
		 		
		double rm = Macrophage.getMres(context).size(); 
		double n = Neutrophil.getNeutrophils(context).size();
		double na = Neutrophil.getApoptosed(context).size();
		double m1ae = Macrophage.getM1ae(context).size();
		double m1de = Macrophage.getM1de(context).size();
		double m1Mac = ((Macrophage.getM1(context).size() - m1de - m1ae));
		double m2Mac = Macrophage.getM2(context).size();
		  
//		// Growth factor secretion calculations - healthy muscle
//		double dtgfdt = (1 * numActiveFibrob + 1 * m1ae + 2 * m2Mac); // tgf
//		double dtnfdt = (1 * resMac + 2 * n + 2 * m1Mac + 2 * m1ae + 2 * m1de + .2 * inflamSSC); // tnf
//		double digf1dt = (2 * numActiveFibrob + 1 * m1Mac + 1 * m1ae + 1 * m1de + 1 * m2Mac); // igf1
//		double dpdgfdt = (1 * numActiveFibrob); // pdgf
//		double dmmpxdt = (1 * numActiveFibrob + 1 * numActiveSecretingSSC + 1 * numMyofbs);// mmpX
//		double dactiveTGFTempdt = (2 * numMyofbs);// tgf myofibroblast release
//		double decmprotdt = (2 * numActiveFibrob + 1 * numActiveSecretingSSC + 1 * numMyofbs);// ecmprot
//		double dil1dt = (2 * resMac + 1 * n + 1 * m1Mac + 1 * m1ae + 1 * m1de+ 1 * inflamFibroblasts + 1 * numActiveSecretingSSC);// il1
//		double dil8dt = (1 * resMac + 1 * n + 1 * m1Mac + 2 * m1de + 1 * inflamFibroblasts + .2 * inflamSSC);// il8
//		double dcxcl2dt = (1 * resMac);// cxcl2
//		double dcxcl1dt = (1 * resMac);// cxcl1
//		double dccl3dt = (1 * resMac + 1 * n - 1 * na);// ccl3
//		double dccl4dt = (1 * resMac + 1 * n);// ccl4
//		double dil6dt = (1 * n + 3 * m1Mac + 3 * m1de + 1 * numActiveFibrob + 1 * inflamSSC);// il6
//		double dmcpdt = (1.7 * n + 1 * inflamFibroblasts + 1 * inflamSSC);// mcp
//		double difndt = (5 * n + 1 * m1de);// ifn
//		double dlactoferinsdt = (5 * na);// lactoferins
//		double dhgfdt = percentNecrotic * 100 * 5;// hgf - Released from ecm with damage // eliminated effect of apop neutrophils and released with % necrotic at damage
//		double dvegfdt = (1 * na + 1 * m1Mac + m1ae + 1 * m1de + .5 * numActiveSecretingSSC);// vegf
//		double dmmp12dt = (2 * m1Mac + 2 * m1ae + 2 * m1de + 1 * numActiveSecretingSSC);// mmp12
//		double dgcsfdt = (1 * m1Mac + 1 * m1ae + 1 * m1de);// gcsf
//		double dil10dt = (1 * m1Mac + 0.1 * m1ae + 0.5 * m1de + 1.2 * m2Mac);// il10
//		double dlipoxinsdt = (1 * m1ae + 1 * m1de);// lipoxins
//		double dresolvinsdt = (3 * m1ae + 2 * m1de);// resolvins
//		double dccl17dt = (1 * m1ae + 2.8 * m2Mac);// ccl17
//		double dccl22dt = (1 * m1ae + 1 * m2Mac + 1 * numActiveSecretingSSC);// ccl22
//		double dcollagen4dt = (1 * m2Mac);// collagen4
//		double dpge2dt = (3 * m2Mac);// pge2
//		double drosdt = (1 * n + 1 * m1de);// ros
//		double dfgfdt = numActiveFibrob; // fgf
//		double dil4dt = m2Mac; // il4		
		
		
		// Growth factor secretion calculations - dystrophic muscle
		//double dtgfdt = (6 * numActiveFibrob + 1 * m1ae + 2 * m2Mac); // tgf
		double dtgfdt = (1 * numActiveFibrob + 1 * m1ae + 2 * m2Mac); // tgf
		double dtnfdt = (1 * resMac + 2 * n + 2 * m1Mac + 2 * m1ae + 2 * m1de + .2 * inflamSSC)*2; // tnf
		double digf1dt = (2 * numActiveFibrob + 1 * m1Mac + 1 * m1ae + 1 * m1de + 1 * m2Mac); // igf1
		double dpdgfdt = (1 * numActiveFibrob); // pdgf
		double dmmpxdt = (0.2 * numActiveFibrob + 1 * numActiveSecretingSSC + 5 * numMyofbs);// mmpX
		double dactiveTGFTempdt = (2 * numMyofbs);// tgf myofibroblast release
		double decmprotdt = (10 * numActiveFibrob + 1 * numActiveSecretingSSC + 1 * numMyofbs);// ecmprot
		double dil1dt = (2 * resMac + 1 * n + 1 * m1Mac + 1 * m1ae + 1 * m1de + 1 * inflamFibroblasts + 1 * numActiveSecretingSSC);// il1
		double dil8dt = (1 * resMac + 1 * n + 1 * m1Mac + 2 * m1de + 1 * inflamFibroblasts + .2 * inflamSSC);// il8
		double dcxcl2dt = (1 * resMac);// cxcl2
		double dcxcl1dt = (1 * resMac);// cxcl1
		double dccl3dt = (1 * resMac + 1 * n - 1 * na);// ccl3
		double dccl4dt = (1 * resMac + 1 * n);// ccl4
		double dil6dt = (1 * n + 3 * m1Mac + 3 * m1de + 1 * numActiveFibrob + 1 * inflamSSC);// il6
		double dmcpdt = (1.7 * n + 1 * inflamFibroblasts + 1 * inflamSSC);// mcp
		double difndt = (5 * n + 1 * m1de + 1 * m1Mac);// ifn
		double dlactoferinsdt = (5 * na);// lactoferins
		double dhgfdt = percentNecrotic * 100 * 5 + na;// hgf - Released from ecm with damage // eliminated effect of apop neutrophils and released with % necrotic at damage
		double dvegfdt = (1 * na + 1 * m1Mac + m1ae + 1 * m1de + .5 * numActiveSecretingSSC);// vegf
		double dmmp12dt = (2 * m1Mac + 2 * m1ae + 2 * m1de + 1 * numActiveSecretingSSC);// mmp12
		double dgcsfdt = (1 * m1Mac + 1 * m1ae + 1 * m1de);// gcsf
		double dil10dt = (1 * m1Mac + 0.1 * m1ae + 0.5 * m1de + 1.2 * m2Mac);// il10
		double dlipoxinsdt = (1 * m1ae + 1 * m1de);// lipoxins
		double dresolvinsdt = (3 * m1ae + 2 * m1de);// resolvins
		double dccl17dt = (1 * m1ae + 2.8 * m2Mac);// ccl17
		double dccl22dt = (1 * m1ae + 1 * m2Mac + 1 * numActiveSecretingSSC);// ccl22
		double dcollagen4dt = (1 * m2Mac);// collagen4
		double dpge2dt = (3 * m2Mac);// pge2
		double drosdt = (1 * n + 1 * m1de);// ros
		double dfgfdt = numActiveFibrob; // fgf
		double dil4dt = m2Mac; // il4				

		// GROWTH FACTOR SOLVER
		// add new growth factors and include half-life
		growthFactorsTemp[0] = (growthFactors[0] + dtgfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[1] = (growthFactors[1] + dtnfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[2] = (growthFactors[2] + digf1dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[3] = (growthFactors[3] + dpdgfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[4] = (growthFactors[4] + dmmpxdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[5] = (growthFactors[5] + dactiveTGFTempdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[6] = (growthFactors[6] + decmprotdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[7] = (growthFactors[7] + dil1dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[8] = (growthFactors[8] + dil8dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[9] = (growthFactors[9] + dcxcl2dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[10] = (growthFactors[10] + dcxcl1dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[11] = (growthFactors[11] + dccl3dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[12] = (growthFactors[12] + dccl4dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[13] = (growthFactors[13] + dil6dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[14] = (growthFactors[14] + dmcpdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[15] = (growthFactors[15] + difndt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[16] = (growthFactors[16] + dlactoferinsdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[17] = (growthFactors[17] + dhgfdt) * Math.pow(.5, timestep / 8.); // slower half-life to
																							// account for ecm
																							// breakdown/release
		growthFactorsTemp[18] = (growthFactors[18] + dvegfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[19] = (growthFactors[19] + dmmp12dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[20] = (growthFactors[20] + dgcsfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[21] = (growthFactors[21] + dil10dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[22] = (growthFactors[22] + dlipoxinsdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[23] = (growthFactors[23] + dresolvinsdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[24] = (growthFactors[24] + dccl17dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[25] = (growthFactors[25] + dccl22dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[26] = (growthFactors[26] + dcollagen4dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[27] = (growthFactors[27] + dpge2dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[28] = drosdt; // does not build up- flushes every step
		growthFactorsTemp[29] = (growthFactors[29] + dfgfdt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[30] = (growthFactors[30] + dil4dt) * Math.pow(.5, timestep / 5.);
		growthFactorsTemp[31] = growthFactors[31];
		growthFactors = growthFactorsTemp;
		// if growth factor below threshold == 0
		if (growthFactors[0] < .0001) {
			growthFactors[0] = 0.;
		}
		if (growthFactors[1] < .0001) {
			growthFactors[1] = 0.;
		}
		
//		// Full TNF block
		// growthFactors[1] = 0.;
//		// 
		if (growthFactors[2] < .0001) {
			growthFactors[2] = 0.;
		}
		if (growthFactors[3] < .0001) {
			growthFactors[3] = 0.;
		}
		if (growthFactors[4] < .0001) {
			growthFactors[4] = 0.;
		}
		if (growthFactors[5] < .0001) {
			growthFactors[5] = 0.;
		}
		if (growthFactors[6] < .0001) {
			growthFactors[6] = 0.;
		}
		if (growthFactors[7] < .0001) {
			growthFactors[7] = 0.;
		}
		if (growthFactors[8] < .0001) {
			growthFactors[8] = 0.;
		}
		if (growthFactors[9] < .0001) {
			growthFactors[9] = 0.;
		}
		if (growthFactors[10] < .0001) {
			growthFactors[10] = 0.;
		}
		if (growthFactors[11] < .0001) {
			growthFactors[11] = 0.;
		}
		if (growthFactors[12] < .0001) {
			growthFactors[12] = 0.;
		}
		if (growthFactors[13] < .0001) {
			growthFactors[13] = 0.;
		}
		if (growthFactors[14] < .0001) {
			growthFactors[14] = 0.;
		}
		if (growthFactors[15] < .0001) {
			growthFactors[15] = 0.;
		}
		if (growthFactors[16] < .0001) {
			growthFactors[16] = 0.;
		}
		if (growthFactors[17] < .0001) {
			growthFactors[17] = 0.;
		}
		if (growthFactors[18] < .0001) {
			growthFactors[18] = 0.;
		}
		if (growthFactors[19] < .0001) {
			growthFactors[19] = 0.;
		}
		if (growthFactors[20] < .0001) {
			growthFactors[20] = 0.;
		}
		if (growthFactors[21] < .0001) {
			growthFactors[21] = 0.;
		}
		if (growthFactors[22] < .0001) {
			growthFactors[22] = 0.;
		}
		if (growthFactors[23] < .0001) {
			growthFactors[23] = 0.;
		}
		if (growthFactors[24] < .0001) {
			growthFactors[24] = 0.;
		}
		if (growthFactors[25] < .0001) {
			growthFactors[25] = 0.;
		}
		if (growthFactors[26] < .0001) {
			growthFactors[26] = 0.;
		}
		if (growthFactors[27] < .0001) {
			growthFactors[27] = 0.;
		}
		if (growthFactors[28] < .0001) {
			growthFactors[28] = 0.;
		}
		if (growthFactors[29] < .0001) {
			growthFactors[29] = 0.;
		}
		if (growthFactors[30] < .0001) {
			growthFactors[30] = 0.;
		}
		return growthFactors;
	}

	// Calculate active-TGFbeta
	public static void initializeActiveTgf() {
		for (int i = 0; i < activeDelay; i++) {
			activeTgf[i] = 0 + mdxTGF;
		}
	}

	public static void setActiveTgf(int tick) {
		// at each time step, store the amount of active tgf
		activeTgf[tick + activeDelay] = growthFactors[0] + mdxTGF;
		// activeTgf[tick + activeDelay] = 0; // TEST: block active tgf
	}

	public static double getActiveTgf(int tick) {
		return activeTgf[tick] + growthFactors[5];
		// return 0; // TEST: block active tgf
	}

	// get growth factors:
	public double getTgf() {
		return growthFactors[0];
	}

	public double getTnf() {
		return growthFactors[1] + Fiber.mdxTnf;
	}

	public double getIgf1() {
		return growthFactors[2];
	}

	public double getPdgf() {
		return growthFactors[3];
	}

	public double getMmpx() {
		return growthFactors[4];
	}

	public double getCollagenX() {
		return growthFactors[5];
	}

	public double getEcmprot() {
		return growthFactors[6];
	}

	public double getIl1() {
		return growthFactors[7];
	}

	public double getIl8() {
		return growthFactors[8];
	}

	public double getCxcl2() {
		return growthFactors[9];
	}

	public double getCxcl1() {
		return growthFactors[10];
	}

	public double getCcl3() {
		return growthFactors[11];
	}

	public double getCcl4() {
		return growthFactors[12];
	}

	public double getIl6() {
		return growthFactors[13];
	}

	public double getMcp() {
		return growthFactors[14];
	}

	public double getIfn() {
		return growthFactors[15] + Fiber.mdxIfn;
	}

	public double getLactoferins() {
		return growthFactors[16];
	}

	public double getHgf() {
		return growthFactors[17];
	}

	public double getVegf() {
		return growthFactors[18];
	}

	public double getMmp12() {
		return growthFactors[19];
	}

	public double getGcsf() {
		return growthFactors[20];
	}

	public double getIl10() {
		return growthFactors[21];
	}

	public double getLipoxins() {
		return growthFactors[22];
	}

	public double getResolvins() {
		return growthFactors[23];
	}

	public double getCcl17() {
		return growthFactors[24];
	}

	public double getCcl22() {
		return growthFactors[25];
	}

	public double getCollagen4() {
		return growthFactors[26];
	}

	public double getPge2() {
		return growthFactors[27];
	}

	public double getRos() {
		return growthFactors[28];
	}

	public double getFgf() {
		return growthFactors[29];
	}
	
	public double getIl4() {
		return growthFactors[30];
	}
	
	public double getInflamWeight() {
		return inflamWeight;
	}
	
	public double getTotalMCP() {
		return growthFactors[31];
	}
	
	public double mcpRatio() {
		return growthFactors[31]/growthFactors[14];
	}
	

}
