/**
 * 
 */
package muscle;

import java.util.ArrayList;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;
import repast.simphony.valueLayer.GridValueLayer;

/**
 * @author Kelley Virgilio This class includes the inflammatory cell
 *         calculations for cells NOT spatially located
 * 
 */
public class InflamCell {

	private static Grid<Object> grid;
	public static ContinuousSpace<Object> space;
	private static BufferedGridValueLayer mcpSpatial;
	public static double[] inflamCellsIter; // Use inflamCellsIter to be able to solve ODE at each step
	static final int numInflamCells = 10; // Number of inflammatory cells tracked in simulations
	public static double[] inflamCells = new double[numInflamCells]; // Array of inflammatory cell counts = [RM N Na M1
																		// M1ae M1de M2]
	public static double rmBasal;
	//public static int tick = 0; // keep track of tick count for the active tgf

	// DISEASE STATE PARAMETERS
	public static final int chronicDamage = 1; // if 1 = chronic damage, otherwise single level
	public static double mdxChronicInflam = 1; // 1 at healthy control--> increases number of resident macs

	public InflamCell(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, ContinuousSpace<Object> space) {
		this.grid = grid;
		this.space = space;
		this.mcpSpatial = mcpSpatial;
	}

	// GLOBAL STEPS FOR EACH TIME POINT: SINGLE UPDATE PARAMETER CHANGES
	@ScheduledMethod(start = 2, interval = 1, priority = 1) // only 1 inflammatory cell agent
	public void scheduler() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		// UPDATE INFLAMMATORY CELLS AND GROWTH FACTORS
		double[] inflamCellsIter = getInflamCellsIter(); //
		double[] inflamCells = InflamCell.cellCountSolver(inflamCellsIter, context, Fiber.origFiberNumber); // AT EACH
																											// STEP
		double[] growthFactors = GrowthFactors.growthFactorSolver(inflamCells, context); // solve for growth factor
		setInflamCellsIter(inflamCells); // reset inflamCellsIter for next generation
		GrowthFactors.setActiveTgf(getTick());
		int totalFiberNumber = Fiber.getTotalFiberNumber(context);
		Necrosis.necrosisBehaviors(context, grid, inflamCells, totalFiberNumber, growthFactors);
		setTick();
	}

	// CHRONIC DAMAGE-- Repetitive injury, followed by single damage
	//@ScheduledMethod(start = 1.9, interval = 1)
//	public void chronicDamageScheduleDouble(){
//		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
//		if (tick < 240) {
//			if (chronicDamage == 1 && tick == 0 || chronicDamage == 1 && tick % 24 == 0){ // determines multiples of 24 by finding remainder-- micro damage at 2 weeks then stop
//				Necrosis.chronicDamage(context, grid, 2);
//			}
//		}
//		else if (chronicDamage == 1 && tick == 240){
//			Necrosis.necrosisInitialize(context, grid, 20);
//		}
//	}
	
	// DOUBLE INJURY - Once at start and again at 28 days
			//@ScheduledMethod(start = 1.9, interval = 1)
			/*public void chronicDamageScheduleDouble() {
				Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
				if (Fiber.tick == 0) {
					if (chronicDamage == 1 && Fiber.tick == 0 || chronicDamage == 1 && Fiber.tick == 672) { // determines multiples of 24 by finding remainder-- micro damage
						Necrosis.chronicDamage(context, grid, 10);
						// Check every time there is damage:
						// Define the amount of necrosis/fiber (similar to originalFiberNecrosis for the single damage simulations:
						double[] chronicFiberNecrosisTemp = new double[Fiber.origFiberNumber];
						for (int i = 1; i < Fiber.origFiberNumber + 1; i++) { // go through each fiber and change the border to red
							// get a random fiber in order to call getFiberBorder
							List<Object> elemsInFiber = Fiber.getElemInFiber(i, context); // get elems within this fiber
							if (elemsInFiber.size() > 0) {
								Object randomFiber = elemsInFiber.get(0);// choose one fiber to call getFiberBorder
								((Fiber) randomFiber).getFiberBorder(i, context); // get all the borders and set to 1
								for (Object elems : elemsInFiber) {
									if (((Fiber) elems).getDamaged() != 0) {
										chronicFiberNecrosisTemp[i - 1] = 1; // if any of the fibers are marked as damaged- end and
																				// go to the next fiber and check
									}
								}
							}
							Fiber.chronicFiberNecrosis = chronicFiberNecrosisTemp;
						}
					}
				} 
			}
	*/
	// CHRONIC DAMAGE-- Repetitive daily injury - variable with age
		@ScheduledMethod(start = 1.9, interval = 1)
		public void chronicDamageVariableSchedule() {
			Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
			if (Fiber.tick <= 336) {
				if (chronicDamage == 1 && Fiber.tick == 0 || chronicDamage == 1 && Fiber.tick % 24 == 0) { // determines multiples of 24 by finding remainder-- micro damage
					Necrosis.chronicDamage(context, grid, 2);
					// Check every time there is damage:
					// Define the amount of necrosis/fiber (similar to originalFiberNecrosis for the single damage simulations:
					double[] chronicFiberNecrosisTemp = new double[Fiber.origFiberNumber];
					for (int i = 1; i < Fiber.origFiberNumber + 1; i++) { // go through each fiber and change the border to red
						// get a random fiber in order to call getFiberBorder
						List<Object> elemsInFiber = Fiber.getElemInFiber(i, context); // get elems within this fiber
						if (elemsInFiber.size() > 0) {
							Object randomFiber = elemsInFiber.get(0);// choose one fiber to call getFiberBorder
							((Fiber) randomFiber).getFiberBorder(i, context); // get all the borders and set to 1
							for (Object elems : elemsInFiber) {
								if (((Fiber) elems).getDamaged() != 0) {
									chronicFiberNecrosisTemp[i - 1] = 1; // if any of the fibers are marked as damaged- end and
																			// go to the next fiber and check
								}
							}
						}
						Fiber.chronicFiberNecrosis = chronicFiberNecrosisTemp;
					}
				}
			} else if (chronicDamage == 1 && Fiber.tick % 24 == 0) {
				// fading damage as mouse ages (Evans et al. 2010)
				double damage = 1; 
				if (Fiber.tick > 2496 && Fiber.tick < 3984) { //3-6 months less damage
					damage = 0.5; 
				} else if (Fiber.tick >= 3984) { //less damage after 6 months for rest of life 
					damage = 0.25;
				}
				Necrosis.chronicDamage(context, grid, damage);
				// Check every time there is damage:
				// Define the amount of necrosis/fiber (similar to originalFiberNecrosis for the single damage simulations:
				double[] chronicFiberNecrosisTemp = new double[Fiber.origFiberNumber];
				for (int i = 1; i < Fiber.origFiberNumber + 1; i++) { // go through each fiber and change the border to red
					// get a random fiber in order to call getFiberBorder
					List<Object> elemsInFiber = Fiber.getElemInFiber(i, context); // get elems within this fiber
					if (elemsInFiber.size() > 0) {
						Object randomFiber = elemsInFiber.get(0);// choose one fiber to call getFiberBorder
						((Fiber) randomFiber).getFiberBorder(i, context); // get all the borders and set to 1
						for (Object elems : elemsInFiber) {
							if (((Fiber) elems).getDamaged() != 0) {
								chronicFiberNecrosisTemp[i - 1] = 1; // if any of the fibers are marked as damaged- end and
																		// go to the next fiber and check
							}
						}
					}
					Fiber.chronicFiberNecrosis = chronicFiberNecrosisTemp;
				}
			}
		}
	
	// CHRONIC DAMAGE-- Repetitive injury, same amount every day
	//@ScheduledMethod(start = 1.9, interval = 1)
	/*public void chronicDamageSchedule() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		if (chronicDamage == 1 && Fiber.tick == 0 || chronicDamage == 1 && Fiber.tick % 24 == 0) { // determines multiples of 24 by finding remainder-- micro damage
			Necrosis.chronicDamage(context, grid, 0.5); // change injury amount here (percentage)
			// Check every time there is damage:
			// Define the amount of necrosis/fiber (similar to originalFiberNecrosis for the single damage simulations:
			double[] chronicFiberNecrosisTemp = new double[Fiber.origFiberNumber];
			for (int i = 1; i < Fiber.origFiberNumber + 1; i++) { // go through each fiber and change the border to red
				// get a random fiber in order to call getFiberBorder
				List<Object> elemsInFiber = Fiber.getElemInFiber(i, context); // get elems within this fiber
				if (elemsInFiber.size() > 0) {
					Object randomFiber = elemsInFiber.get(0);// choose one fiber to call getFiberBorder
					((Fiber) randomFiber).getFiberBorder(i, context); // get all the borders and set to 1
					for (Object elems : elemsInFiber) {
						if (((Fiber) elems).getDamaged() != 0) {
							chronicFiberNecrosisTemp[i - 1] = 1; // if any of the fibers are marked as damaged- end and
																	// go to the next fiber and check
						}
					}
				}
				Fiber.chronicFiberNecrosis = chronicFiberNecrosisTemp;
			}
		}

	}
	
//	// CHRONIC DAMAGE-- Repetitive injury
	//@ScheduledMethod(start = 1.9, interval = 1)
	public void chronicDamageSchedule2(){ //every 120 hours
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		if (chronicDamage == 1 && Fiber.tick == 0 || chronicDamage == 1 && Fiber.tick % 120 == 30){ // determines multiples of 24 by finding remainder-- micro damage
			Necrosis.chronicDamage(context, grid, 1);
			// Check every time there is damage:
			// Define the amount of necrosis/fiber (similar to originalFiberNecrosis for the single damage simulations:
			double[] chronicFiberNecrosisTemp = new double[Fiber.origFiberNumber];
			for (int i = 1; i < Fiber.origFiberNumber + 1; i++){ // go through each fiber and change the border to red
		    	// get a random fiber in order to call getFiberBorder
		    	List<Object> elemsInFiber = Fiber.getElemInFiber(i, context); // get elems within this fiber
		    	if (elemsInFiber.size() > 0){
		    		Object randomFiber = elemsInFiber.get(0);// choose one fiber to call getFiberBorder
			    	((Fiber) randomFiber).getFiberBorder(i, context); // get all the borders and set to 1
					for (Object elems : elemsInFiber){
						if (((Fiber) elems).getDamaged() != 0){
							chronicFiberNecrosisTemp[i-1] = 1; // if any of the fibers are marked as damaged- end and go to the next fiber and check
						}
					}
				}
			    Fiber.chronicFiberNecrosis = chronicFiberNecrosisTemp;
			}
		}
	}
	
*/
	public static void setTick() {
		Fiber.tick = Fiber.tick + 1;
	}

	public static int getTick() {
		return Fiber.tick;
	}

	static final int timestep = 1; // temp value for the timestep = 1 hour

	public static double[] initialize(int origFiberNumber, double mdxChronicInflam, Context<Object> context) {
		// initalize inflammatory cells and return the counts
		Parameters params = RunEnvironment.getInstance().getParameters(); // RunEnvironment --> provides access to the
																			// environment in which a particular model
																			// runs
		
		inflamCells[0] = 0;
		inflamCells[1] = (Integer) params.getValue("n_initial");
		inflamCells[2] = (Integer) params.getValue("na_initial");
		inflamCells[3] = (Integer) params.getValue("m1_initial");
		inflamCells[4] = (Integer) params.getValue("m1ae_initial");
		inflamCells[5] = (Integer) params.getValue("m1de_initial");
		inflamCells[6] = (Integer) params.getValue("m2_initial");
		return inflamCells;
	}
	

	public static double[] cellCountSolver(double[] inflamCellsIter, Context<Object> context, int origFiberNumber) {
		// Solve for the number of inflammatory cells
		// Take in inflamCells from the previous iteration

		// GET FIBROBLAST, SSC COUNTS AND AMOUNT OF NECROSIS
		double[] inflamCellsTemp = new double[numInflamCells];
		int numFactive = Fibroblast.getActiveFibroblasts(context).size(); // number of active fibroblasts
		int numActiveSecretingSSC = SSC.getNumActSecretingSSCs(context); // number of actively secreting sscs
		double percentNecrotic = Necrosis.getPercentNecrotic(context); // Get percent necrotic of muscle
		double necroticCheck = 0.;
		if (percentNecrotic > 0.001) {
			necroticCheck = 1; // else necroticCheck = 0 which means there is no more necrosis-- instead of
								// percent Necrotic
		}

		// SOLVE FOR INFLAMMATORY CELL COUNTS

		// Resident macrophages
		//double dRMdt = (.1 * rmBasal * .3 - .1 * inflamCellsIter[0]) * 1.8;
		double dRMdt = Macrophage.getMres(context).size(); 
	
		double na = Neutrophil.getApoptosed(context).size();
		double m1 = Macrophage.getM1(context).size();
		double m1ae = Macrophage.getM1ae(context).size();
		double m1de = Macrophage.getM1de(context).size();
		double m2 = Macrophage.getM2(context).size();
		
		// Neutrophils
		// Neutrophils need to arrive at all damage due to broken blood vessels:
		/*double dNdt = ((55.3 * inflamCellsIter[0] * percentNecrotic * 1.3 + .01 * numFactive
				+ .6048 * inflamCellsIter[1] - .72 * inflamCellsIter[2] + (-4.889) * inflamCellsIter[3] * .5
				- .59 * inflamCellsIter[4] - .565 * inflamCellsIter[5] + (-1.34) * inflamCellsIter[6])
				+ .07 * numActiveSecretingSSC) * (1. / 24.) * .8
				+ (-.2 * inflamCellsIter[1] * percentNecrotic * 1.5 * 7 - .1 * inflamCellsIter[1]) * (1. / 24.) * 1.44;
				
		if (Necrosis.getInitialBurstNecrotic(context) > 0 && dNdt < 0) {
			// if there was recent damage but dNdt == 0 --> for instance at chronic damage--
			// then add
			dNdt = dNdt + Necrosis.getInitialBurstNecrotic(context) * 55.3 * inflamCellsIter[0] * 1.3;
		}
		*/
		double dNdt = Neutrophil.getNeutrophils(context).size();

		// Apoptotic neutrophils
		//double dNadt = (.2 * inflamCellsIter[1] * percentNecrotic * 1.5 * 7
			//	- .1 * inflamCellsIter[3] * inflamCellsIter[2] - .1 * inflamCellsIter[2]) * (1. / 24.) * 1.8 * .8;
		double dNadt = Neutrophil.getApoptosed(context).size();

		
		// M1 macrophages
		/*double dM1dt = ((.135 * inflamCellsIter[0] * percentNecrotic - .054 * numFactive
				+ .661 * inflamCellsIter[1] * 1.5 - .054 * inflamCellsIter[2] + (.2349) * inflamCellsIter[3]
				- .1647 * inflamCellsIter[4] + .2889 * inflamCellsIter[5] + (-.3483) * inflamCellsIter[6] * .8
				+ .054 * numActiveSecretingSSC) * (1. / 24.) * .8
				+ (-.1 * inflamCellsIter[3] * inflamCellsIter[2] * .5 - .1 * inflamCellsIter[3] * necroticCheck * 2
						- .1 * inflamCellsIter[3] * Fiber.mdxM1death) * (1. / 24.) * .54)
				* Fiber.mdxM1mult;
				*/
		double dM1dt = Macrophage.getM1(context).size();

		// M1- apoptotic neutrophil eating
		//double dM1aedt = (.1 * inflamCellsIter[3] * inflamCellsIter[2] * .5 - .1 * inflamCellsIter[4] * 2
			//	- .1 * inflamCellsIter[4]) * (1. / 24.) * 1.8 * .8;
		double dM1aedt = Macrophage.getM1ae(context).size();

		// M1- debris eating
		//double dM1dedt = (.1 * inflamCellsIter[3] * necroticCheck * 2 - .1 * inflamCellsIter[5] * 2) * (1. / 24.) * 1.8
			//	* .8;
		double dM1dedt = Macrophage.getM1de(context).size();
		
		// M2 macrophages
		/*double dM2dt = ((.18 * inflamCellsIter[0] * percentNecrotic + -.09 * numFactive - .027 * inflamCellsIter[1]
				- .09 * inflamCellsIter[2] + (.7398) * inflamCellsIter[3] * 2 - .126 * inflamCellsIter[4] * .8
				+ (.3456) * inflamCellsIter[5] * 1.5 + (-.18) * inflamCellsIter[6] + 0.09 * numActiveSecretingSSC)
				* (1. / 24.) * .7
				+ (.1 * inflamCellsIter[4] * 2 - .1 * inflamCellsIter[6] * Fiber.mdxM2death) * (1. / 24.)) * .9
				* Fiber.mdxM2mult;
				*/
		double dM2dt = Macrophage.getM2(context).size();

		/*
		// Solve for new counts
		inflamCellsTemp[0] = inflamCellsIter[0] + dRMdt;
		if (inflamCellsIter[1] + dNdt > origFiberNumber * 2) { // Neutrophil ceiling for probing inflammatory cells
			inflamCellsTemp[1] = inflamCellsIter[1];
		} else {
			inflamCellsTemp[1] = inflamCellsIter[1] + dNdt;
		}
		inflamCellsTemp[2] = inflamCellsIter[2] + dNadt;
		inflamCellsTemp[3] = inflamCellsIter[3] + dM1dt;
		inflamCellsTemp[4] = inflamCellsIter[4] + dM1aedt;
		inflamCellsTemp[5] = inflamCellsIter[5] + dM1dedt;
		inflamCellsTemp[6] = inflamCellsIter[6] + dM2dt;
		*/

		// Solve for new counts
		inflamCellsTemp[0] = dRMdt;
		//System.out.println("Res Macs: ");
		//System.out.println(inflamCellsTemp[0]);
		inflamCellsTemp[1] = dNdt;
		inflamCellsTemp[2] = dNadt;
		inflamCellsTemp[3] = dM1dt;
		inflamCellsTemp[4] = dM1aedt;
		inflamCellsTemp[5] = dM1dedt;
		inflamCellsTemp[6] = dM2dt;
		
		// Calculate fibroblast and ssc additional parameters:
		// FIBROBLAST RECRUITMENT- EOSINOPHIL-IL4
		double dIL4Recruitdt = Necrosis.getInitialBurstNecrotic(context) * 10 * origFiberNumber * Fibroblast.fibroScale * 30 * Fibroblast.murphRepRecruitParam * Fiber.mdxEosinophil * Fiber.tcf4Scale;
		inflamCellsTemp[7] = (inflamCellsIter[7] + dIL4Recruitdt) * Math.pow(.5, timestep / (24.)); // slow decrease to allow necrosis signal to be maintained for a long time
		// FIBROBLAST SSC-DEPENDENT EXPANSION
		double numRecruitSSC = SSC.getActiveSSCs(context).size() - SSC.getDiffSSCs(context).size();// number of active, not differentiated sscs
		double sscFibroRecruit = numRecruitSSC / Fiber.getTotalFiberNumber(context); // weight by the number of fibers
		double dfibrobExpansiondt = 200 * sscFibroRecruit ; // fibroblast recruitment parameter
		inflamCellsTemp[8] = dfibrobExpansiondt;
				//(inflamCellsIter[8] + dfibrobExpansiondt) * Math.pow(.5, timestep / 5.);
		// SSC ACTIVATION
		// Use instead of just hgf, since hgf is dependent on how quickly damage is cleared
		//double dsscActivationdt = Necrosis.getInitialBurstNecrotic(context) * 10 * origFiberNumber * 15;
		//inflamCellsTemp[9] = (inflamCellsIter[9] + dsscActivationdt) * Math.pow(.5, timestep / 24.); // slow decrease to allow necrosis signal to be maintained for a long time
//		double dsscActivationdt = GrowthFactors.growthFactors[17] + GrowthFactors.growthFactors[17] + GrowthFactors.growthFactors[2]; // 2*hgf + fgf + igf
//		inflamCellsTemp[9] = dsscActivationdt ; // slow decrease to allow necrosis signal to be maintained for a long time
//		
		//System.out.println("SSC: " + inflamCellsTemp[9]);
		//System.out.println("hgf: " + GrowthFactors.growthFactors[17]);
		//System.out.println("fgf: " + GrowthFactors.growthFactors[29]);


		// Check to make sure all numbers are not less than zero/lower threshold
		if (inflamCellsTemp[0] < .001) {
			inflamCellsTemp[0] = 0.;
		}
		if (inflamCellsTemp[1] < .001) {
			inflamCellsTemp[1] = 0.;
		}
		if (inflamCellsTemp[2] < .001) {
			inflamCellsTemp[2] = 0.;
		}
		if (inflamCellsTemp[3] < .001) {
			inflamCellsTemp[3] = 0.;
		}
		if (inflamCellsTemp[4] < .001) {
			inflamCellsTemp[4] = 0.;
		}
		if (inflamCellsTemp[5] < .001) {
			inflamCellsTemp[5] = 0.;
		}
		if (inflamCellsTemp[6] < .001) {
			inflamCellsTemp[6] = 0.;
		}
		if (inflamCellsTemp[7] < .001) {
			inflamCellsTemp[7] = 0.;
		}
		if (inflamCellsTemp[8] < .001) {
			inflamCellsTemp[8] = 0.;
		}
		if (inflamCellsTemp[9] < .001) {
			inflamCellsTemp[9] = 0.;
		}

		// DISEASE STATE PARAMETERS
//		// macrophage knockdown- Shen 2006
//		// TURN OFF AUTO NECROSIS REMOVAL SO THAT IT DOESN'T DISAPPEAR
//		inflamCellsTemp[3] = inflamCellsTemp[3]*.9;
//		inflamCellsTemp[4] = inflamCellsTemp[3]*.9;
//		inflamCellsTemp[5] = inflamCellsTemp[3]*.9;
//		if (getTick()> 20 && getTick()< 40){
//			inflamCellsTemp[1] = inflamCellsTemp[1] - .9*inflamCellsTemp[1];
//			inflamCellsTemp[2] = inflamCellsTemp[2] - .9*inflamCellsTemp[2];
//		}
//		if (getTick()> 40){
//			inflamCellsTemp[1] = 0;
//			inflamCellsTemp[2] = 0;
//		}

//		// MACROPHAGE KNOCKDOWN- at start
//		if (getTick() < 120){
//			inflamCellsTemp[3] = 0;
//			inflamCellsTemp[4] = 0;
//			inflamCellsTemp[5] = 0;
//		}

//		//
//		// macrophage increase- M1 and M2
//		// TURN OFF AUTO NECROSIS REMOVAL SO THAT IT DOESN'T DISAPPEAR
//		inflamCellsTemp[3] = inflamCellsTemp[3]*1.01;
//		inflamCellsTemp[4] = inflamCellsTemp[3]*1.01;
//		inflamCellsTemp[5] = inflamCellsTemp[3]*1.01;
//		//
		// macrophage increase- M2 only
		// TURN OFF AUTO NECROSIS REMOVAL SO THAT IT DOESN'T DISAPPEAR
//		inflamCellsTemp[6] = inflamCellsTemp[6]*1.01;
		//
		inflamCells = inflamCellsTemp;
		return inflamCells;

	}

	public static double[] getInflamCellsIter() {
		return inflamCellsIter;
	}

	public static void setInflamCellsIter(double[] inflamCells) {
		inflamCellsIter = inflamCells;
	}

	public double getRM() {
		return inflamCells[0];
	}

	public double getN() {
		return inflamCells[1];
	}

	public double getNa() {
		return inflamCells[2];
	}

	public double getM1() {
		return inflamCells[3];
	}

	public double getM1Total() {
		return inflamCells[3] + inflamCells[4] + inflamCells[5];
	}

	public double getM1ae() {
		return inflamCells[4];
	}

	public double getM1de() {
		return inflamCells[5];
	}

	public double getM2() {
		return inflamCells[6];
	}

	public static double getFibrobRecruit() {
		return inflamCells[7];
	}

	// Counts and plotting parameters called from inflamCell:
	public double getFibroblastNumber4Plot() { // added to use aggregate counts (only 1 inflam cell, so only called one
												// time)
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		List<Object> fibroblasts = Fibroblast.getFibroblasts(context);// Get total number of fibroblasts
		int numFibroblasts = fibroblasts.size();
		return numFibroblasts;
	}
	
	public double getFAPNumber4Plot() { // added to use aggregate counts (only 1 inflam cell, so only called one
		// time)
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		List<Object> faps = Fibroblast.getFAPs(context);// Get total number of fibroblasts
		int numFAPs = faps.size();
		return numFAPs;
	}

	public double getMyofibroblastNumber4Plot() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		List<Object> myofibroblasts = Fibroblast.getMyofibroblasts(context);// Get total number of fibroblasts
		int numMyofibroblasts = myofibroblasts.size();
		return numMyofibroblasts;
	}

	public double getPax7Number4Plot() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		return (SSC.getSSCs(context).size() - SSC.getDiffSSCs(context).size());
	}

	public double getSenescentSSC4Plot() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		return (SSC.getSenescentSSCs(context));
	}

	public double getMyoblastNumber4Plot() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		return (SSC.getMyoblastSSCs(context).size());
	}

	public double getDiffNumber4Plot() {
		Context context = ContextUtils.getContext(this); // get the context of the inflammatory cell
		return SSC.getDiffSSCs(context).size();
	}

	public double getActiveFibroblast4Plot() {
		Context context = ContextUtils.getContext(this);
		return Fibroblast.getActiveFibroblasts(context).size();
	}

	public int getTotalSize() {
		Context context = ContextUtils.getContext(this);
		return Fiber.getFiberElems(context).size() + ECM.getECM(context).size() + Necrosis.getNecrosis(context).size();
	}

	public double getTotalCollagen() {
		Context context = ContextUtils.getContext(this);
		return ECM.getTotalCollagenAmt(context);
	}

	public double getECMTotal() {
		Context context = ContextUtils.getContext(this);
		return ECM.getECM(context).size();
	}

	public double getFiberTotal() {
		Context context = ContextUtils.getContext(this);
		return Fiber.getFiberElems(context).size();
	}

	public double getNecrosisTotal() {
		Context context = ContextUtils.getContext(this);
		return Necrosis.getNecrosis(context).size();
	}

	public double getCollagenDensity() {
		Context context = ContextUtils.getContext(this);
		ECM.collagenDensity = ECM.getTotalCollagenAmt(context) / ECM.getECM(context).size();
		return ECM.getTotalCollagenAmt(context) / ECM.getECM(context).size();
	}
}
