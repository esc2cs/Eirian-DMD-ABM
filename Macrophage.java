package muscle;

import java.util.ArrayList;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;

public class Macrophage {

	private static ContinuousSpace<Object> space;
	private static Grid<Object> grid;
	private static BufferedGridValueLayer mcpSpatial;
	// private static GridValueLayer mcpSpatial;
	private int phenotype; // 0 = M0, 1 = M1, 2 = M2, 3 = resident
	private int age;
	private int phagocytosis;
	private static double[] growthFactors;
	private Object buddy;
	private boolean linked;
	private int stay;
	private static int m1num; // number of M1s
	private static int m2num; // number of M2s
	private int resMnum; // number of resMs
	private int numPhagocytosed; //temp to test fap phagocytosis efficiency

	public Macrophage(BufferedGridValueLayer mcpSpatial, ContinuousSpace<Object> space, Grid<Object> grid,
			int phenotype, int age, int phagocytosis, boolean linked, Object buddy) {
		this.space = space;
		this.grid = grid;
		this.phenotype = phenotype;
		this.age = age;
		this.phagocytosis = phagocytosis;
		this.linked = linked;
		this.buddy = buddy;
		this.mcpSpatial = mcpSpatial;
		this.numPhagocytosed = 0; //temp to test fap phagocytosis efficiency
	}

	// MACROPHAGE RECRUITMENT
	//@ScheduledMethod(start = 2, interval = 1, pick = 1)
	public static void macRecruitment(Context<Object> context) {

		//Context context = ContextUtils.getContext(this);
		if (Necrosis.getInitialBurstNecrotic(context) > 0) {
			double rmBasal = Math.floor(Fiber.getTotalFiberNumber(context) / 3.7); // defines number of resident macrophages based on number of fibers
			if (rmBasal == 0) {
				rmBasal = 1;
			}
			
			for (int i = 0; i < rmBasal; i++) {
				context.add(new Macrophage(mcpSpatial, space, grid, 3, 0, 0, false, null));
			}
		}
		
		
		double[] growthFactors = GrowthFactors.getGrowthFactors();
		// Recruiting growth factors
		double ccl4 = growthFactors[12];
		double ccl17 = growthFactors[24];
		double ccl22 = growthFactors[25];
		double ccl13 = growthFactors[11];
		double il6 = growthFactors[13];
		double mcp = growthFactors[14];
		// Deterring growth factors
		double pge2 = growthFactors[27];
		double lipoxins = growthFactors[22];
		double tgf = growthFactors[0];
		double ros = growthFactors[28];
		// M1 growth factors
		double ifn = growthFactors[15];
		double tnf = growthFactors[1];
		// M2 growth factors
		double il10 = growthFactors[21];
		double il4 = growthFactors[30];

		double recruit = ccl4 + ccl17 + ccl22 + ccl13 + il6 + mcp; // Azurocidin + LL37 + CCL2 + CCL4 + CCL17 + CCL22 + CCL3 + CCL6 + sum [IL6] of patches + CX3CL1 + sum [MCP] of patches
		double deter = pge2 + lipoxins + tgf; // PGE2 + Lipoxins + NO + sum [TGF] of patches
		double macProb = 0;
		double differential = recruit - deter;

		if (differential > 0) { // threshold to recruit macrophage
			macProb = RandomHelper.nextIntFromTo(0, (int) (15 - 8 + 2 * 6) / 5);/// 10) ; // recruit Macrophages if
																				/// recruitment factors > dettering
		}

		// Determine the likelihood of the recruited macrophage to be m1 or M2
		//double m1Chance = ifn + tnf / 2 - il10 - tgf; // IFN + sum [TNF] of patches / 2 - sum [IL10] of patches - sum[TGF] of patches
		double m1Chance = ifn + tnf;// - il10 - tgf; // IFN + sum [TNF] of patches / 2 - sum [IL10] of patches - sum[TGF] of patches
		if (m1Chance < 0) {
			m1Chance = 0.01;
		}

		//double m2Chance = (il10 * 2 + il4 - ifn); // (sum [IL10] of patches * 2) + IL4 + IL13 - IFN
		double m2Chance = (il10 * 2 + il4 - ifn); // (sum [IL10] of patches * 2) + IL4 + IL13 - IFN
		if (m2Chance < 0) {
			m2Chance = 0.01;
		}

		m1num = getM1(context).size(); // get number of M1s
		m2num = getM2(context).size(); // get number of M2s
		if (m1num > 1000 || m2num > 1000) {
			macProb = 0; // don't recruit Macs if M1 OR M2 > 100 (NetLogo: 1000 instead of 100)
		}
		// Add in macrophages
		while (macProb > 0) {
			double chance = RandomHelper.nextDoubleFromTo(0, m1Chance + m2Chance);
			List<Object> ecms = ECM.getECM(context);
			List<Object> necrosis = Necrosis.getNecrosis(context);
			List<Necrosis> edgeNecrosis = Necrosis.getNecrosisAtEdge(context, grid);
//			List<Object> patch = new ArrayList<Object>();
//			for (Object obj : ecms) {
//				if (obj instanceof ECM && ((ECM) obj).getCollagen() > 0) {
//					patch.add(obj);
//				}
//			}
//			for (Object obj : necrosis) {
//				if (obj instanceof Necrosis) {
//					patch.add(obj);
//				}
//			}

			if (chance < m1Chance) {
				if (edgeNecrosis.size() > 0) {
					int index = RandomHelper.nextIntFromTo(0, edgeNecrosis.size() - 1);
					Object spawnpt = edgeNecrosis.get(index);
					GridPoint gridpt = grid.getLocation(spawnpt);
					// NdPoint spacept = space.getLocation(spawnpt);
					Macrophage m1mac = new Macrophage(mcpSpatial, space, grid, 1, 0, 0, false, null);
					context.add(m1mac);
					grid.moveTo(m1mac, gridpt.getX(), gridpt.getY());
					// space.moveTo(m1mac, spacept.getX(), spacept.getY());
				} else if (necrosis.size() > 0){
					int index = RandomHelper.nextIntFromTo(0, necrosis.size() - 1);
					Object spawnpt = necrosis.get(index);
					GridPoint gridpt = grid.getLocation(spawnpt);
					// NdPoint spacept = space.getLocation(spawnpt);
					Macrophage m1mac = new Macrophage(mcpSpatial, space, grid, 1, 0, 0, false, null); // phenotype = 1; age = 0; phagocytosis = 0; linked = false; buddy = null
					context.add(m1mac);
					grid.moveTo(m1mac, gridpt.getX(), gridpt.getY());
					// space.moveTo(m1mac, spacept.getX(), spacept.getY());
				} else {
					int index = RandomHelper.nextIntFromTo(0, ecms.size() - 1);
					Object spawnpt = ecms.get(index);
					GridPoint gridpt = grid.getLocation(spawnpt);
					// NdPoint spacept = space.getLocation(spawnpt);
					Macrophage m1mac = new Macrophage(mcpSpatial, space, grid, 1, 0, 0, false, null);
					context.add(m1mac);
					grid.moveTo(m1mac, gridpt.getX(), gridpt.getY());
					// space.moveTo(m1mac, spacept.getX(), spacept.getY());
				}
				
			} else {
				int index = RandomHelper.nextIntFromTo(0, ecms.size() - 1);
				Object spawnpt = ecms.get(index);
				GridPoint gridpt = grid.getLocation(spawnpt);
				// NdPoint spacept = space.getLocation(spawnpt);
				Macrophage m2mac = new Macrophage(mcpSpatial, space, grid, 2, 0, 0, false, null); // phenotype = 2; age = 0; phagocytosis = 0; linked = false; buddy = null
				context.add(m2mac);
				grid.moveTo(m2mac, gridpt.getX(), gridpt.getY());
				// space.moveTo(m2mac, spacept.getX(), spacept.getY());
			}
			macProb = macProb - 1;
		}
	}

	// Macrophage step
	@ScheduledMethod(start = 2, interval = 1)
	public void step() {
		// Sense environment
		Context context = ContextUtils.getContext(this);
		growthFactors = GrowthFactors.getGrowthFactors();

		// Age macrophage
		this.age = this.age + 1;

		// Determine behaviors
		if (this.phenotype == 1) { // m1 macrophage
			if (this.phagocytosis == 0) { // hasn't eating anything
				m1Behave(context);
			} else if (this.phagocytosis > 0) { // has eaten apoptotic neutrophil
				m1ApopEat(context);
			} else if (this.phagocytosis < 0) {
				m1DebrisEat(context);
			}
		} else if (this.phenotype == 2) { // m2 macrophage
			m2Behave(context);
		} else if (this.phenotype == 3) {
			// after neutrophils are recruited, resident macrophages will stop being active, but there should always be a pool of ~5 left
			if (Neutrophil.numNeutrophils > 5 && getMres(context).size() > 5 || this.age > 5 && getMres(context).size() > 5) {
				if (RandomHelper.nextIntFromTo(0, 4) == 1) { // originally, NetLogo: random 5
					macDie();
				}
			}
		}

		// Non-resident macrophages proliferate
		if (this.phenotype != 3) {
			// Determine proliferation
			int random = RandomHelper.nextIntFromTo(0, 43); // originally, NetLogo: random 10 + 5 * 6.7
			m1num = getM1(context).size(); // get number of M1s
			m2num = getM2(context).size(); // get number of M2s
			resMnum = getMres(context).size(); // get number of resMs
			int macCount = m1num + m2num + resMnum; // total macrophages
			if ((random == 1)) {// && (macCount < 100)) { // NetLogo: 1000 instead of 100
				// macProliferate();
			}

			// Determine death
			if (this.age > 10 + 5 * 6.2) {
				if (RandomHelper.nextIntFromTo(0, 4) == 1) { // originally, NetLogo: random 5
					macDie();
				}
			}
		}
	}

	// M1 BEHAVIOR
	public void m1Behave(Context<Object> context) {
		// omit M1 secretions

		// Find neutrophil and necrosis in this grid point
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> necroticCell = new ArrayList<Object>();
		List<Object> neutrophilNgh = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			} else if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) { // neutrophil in same
																								// "cell"?
				neutrophilNgh.add(obj);
			}
		}
		
		GridPoint pt = grid.getLocation(this);
		//get the objects on the same location as the macrophage
		Iterable<Object> iterPt = grid.getObjectsAt(pt.getX(), pt.getY());
		for (Object obj : iterPt) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			} else if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) { // neutrophil in same																				// "cell"?
				neutrophilNgh.add(obj);
			}
		}

		// Check on conditions to determine behavior
		
		if (linked == true) { // Macrophages will move towards its buddy. If there is none, its no longer linked 
			if (buddy != null) {
				GridPoint movePt = grid.getLocation(buddy);
				System.out.println("Buddy: ( " + movePt.getX() + " , " + movePt.getY() + " )");
				grid.moveTo(this, movePt.getX(), movePt.getY());
			} else {
				linked = false;
			}
		} else {

			if (necroticCell.size() > 0) {
				if (neutrophilNgh.size() > 0) {
					// choose randomly to each either necrosis (debris) or neutrophil
					if (RandomHelper.nextIntFromTo(0, 1) == 0) { // 0 = eat debris
						for (int i = 0; i < necroticCell.size(); i++) {
							Object necrosisRandom = necroticCell.get(i);
							Necrosis.eatNecrosis(context, grid, necrosisRandom);
							this.numPhagocytosed++;
						}
						this.phagocytosis = phagocytosis - 1;
						if (this.phagocytosis < -8) {
							if (RandomHelper.nextIntFromTo(0, 2) == 1) {
								macDie();
							}
						}
					} else { // 1 = eat neutrophil
						this.phagocytosis = phagocytosis + 1;
						int index = RandomHelper.nextIntFromTo(0, neutrophilNgh.size() - 1);
						Object randNeutrophil = neutrophilNgh.get(index);
						context.remove(randNeutrophil);
					}
				} else { // eat necrosis (debris)
					this.phagocytosis = phagocytosis - 1;
					for (int i = 0; i < necroticCell.size(); i++) {
						Object necrosisRandom = necroticCell.get(i);
						Necrosis.eatNecrosis(context, grid, necrosisRandom);
						this.numPhagocytosed++;
					}
					if (this.phagocytosis < -8) {
						if (RandomHelper.nextIntFromTo(0, 2) == 1) {
							macDie();
						}
					}
				}
			} else { // if not necrosis at current location
				if (neutrophilNgh.size() > 0) {
					// eat neutrophil
					this.phagocytosis = phagocytosis + 1;
					int index = RandomHelper.nextIntFromTo(0, neutrophilNgh.size() - 1);
					Object randNeutrophil = neutrophilNgh.get(index);
					context.remove(randNeutrophil);
				} else {
					m1Migrate();
					if (RandomHelper.nextIntFromTo(0, 24) == 1) { // NetLogo orig: 5 * 4 + 5; chance M1 will polarize to
																	// M2
						stay = 1;
						phenotype = 2;
					}
				}
			}
		}
	}

	public void m1Migrate() {
		Context context = ContextUtils.getContext(this);
		List<Object> necrosisNeighbor = new ArrayList<Object>();
		List<Object> openNeighbor = new ArrayList<Object>();
		List<Object> highCollNeighbor = new ArrayList<Object>();
		double mcpMax = 0;
		GridPoint mcpMax_location = grid.getLocation(this);

		//while (necrosisNeighbor.size() == 0 || i < 4) {
		MooreQuery<Object> query = new MooreQuery<Object>(grid, this, 3, 3); // get neighbors of macrophage
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		necrosisNeighbor = new ArrayList<Object>();
		openNeighbor = new ArrayList<Object>();
		highCollNeighbor = new ArrayList<Object>();

		int high = 0; // this tells us if the neighbor with max mcp value is on high collagen or not
		mcpMax = 0;

		for (Object neighbor : iter) {
			if (neighbor instanceof Necrosis) { // necrotic neighbors
				necrosisNeighbor.add(neighbor);
			} else if (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() > 1) { // neighbors with collagen
				highCollNeighbor.add(neighbor);
				GridPoint pt = grid.getLocation(neighbor);
				double mcp_temp = mcpSpatial.get(pt.getX(), pt.getY());
				if (mcp_temp > mcpMax) {
					mcpMax = mcp_temp;
					high = 1;
					mcpMax_location = pt;
				}
			} else if (neighbor instanceof ECM) {
				openNeighbor.add(neighbor); // ecm neighbors
				GridPoint pt = grid.getLocation(neighbor);
				double mcp_temp = mcpSpatial.get(pt.getX(), pt.getY());
				if (mcp_temp > mcpMax) {
					mcpMax = mcp_temp;
					high = 0;
					mcpMax_location = pt;
				}
			}
		}
		//	i++;
		//}
		

		// if there is necrosis, go there, otherwise move up mcp gradient
		if (necrosisNeighbor.size() > 0) { // if there is necrosis move there
			int index = RandomHelper.nextIntFromTo(0, necrosisNeighbor.size() - 1);
			Object randomNeighbor = necrosisNeighbor.get(index);
			GridPoint pt2 = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt2.getX(), pt2.getY());
		}
		else if (mcpMax > 0) {
		 grid.moveTo(this, mcpMax_location.getX(), mcpMax_location.getY());
		}
		else if (highCollNeighbor.size() > 0) { // if there is any collagen
			int index = RandomHelper.nextIntFromTo(0, highCollNeighbor.size() - 1);
			Object randomNeighbor = highCollNeighbor.get(index);
			GridPoint pt2 = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt2.getX(), pt2.getY());
		} else if (openNeighbor.size() > 0) { // Otherwise just pick a random direction of ecm go to it
			int index = RandomHelper.nextIntFromTo(0, openNeighbor.size() - 1);
			Object randomNeighbor = openNeighbor.get(index);
			GridPoint pt2 = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt2.getX(), pt2.getY());
		}
		else {
		}
	}

	// M1 EATS APOPTOTIC NEUTROPHILS
	public void m1ApopEat(Context<Object> context) {
		// Context context = ContextUtils.getContext(this);
		if (stay > 0) {
			stay = stay + 1;
			age = age - 1;
			if (stay > RandomHelper.nextIntFromTo(0, 2)) { // originally, NetLogo : random (2 * 1.5)
				stay = 0;
				this.phenotype = 2;
			}
		} else {
			MooreQuery<Object> query = new MooreQuery(grid, this, 2, 2); // get "cells" this macrophage is currently in
			Iterable<Object> iter = query.query(); // query the list of agents
			List<Object> apopNeutrophil = new ArrayList<Object>();
			for (Object obj : iter) {
				if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) { // if queried agent is an
																							// apop. neu., add to list
					apopNeutrophil.add(obj);
				}
			}
			
			GridPoint pt = grid.getLocation(this);
			//get the objects on the same location as the macrophage
			Iterable<Object> iterPt = grid.getObjectsAt(pt.getX(), pt.getY());
			for (Object obj : iterPt) {
				if (obj instanceof Neutrophil && ((Neutrophil) obj).getApoptosed() > 0) { // neutrophil in same																				// "cell"?
					apopNeutrophil.add(obj);
				}
			}
			
			if (apopNeutrophil.size() > 0) {
				this.phagocytosis = phagocytosis + 1;
				int index = RandomHelper.nextIntFromTo(0, apopNeutrophil.size() - 1);
				Object neutrophil = apopNeutrophil.get(index);
				context.remove(neutrophil);
			} else {
				m1Migrate();
			}
			if (age > (1 + 8 * 4.3)) {
				stay = 1;
				age = 0;
			}
		}
	}

	// M1 EATS DEBRIS
	public void m1DebrisEat(Context<Object> context) {
		// Context context = ContextUtils.getContext(this);
		// secrete cytokine
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get "cells" this macrophage is currently in
		Iterable<Object> iter = query.query(); // query the list of agents
		List<Object> necroticCell = new ArrayList<Object>();
		for (Object obj : iter) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			}
		}
		
		GridPoint pt = grid.getLocation(this);
		//get the objects on the same location as the macrophage
		Iterable<Object> iterPt = grid.getObjectsAt(pt.getX(), pt.getY());
		for (Object obj : iterPt) {
			if (obj instanceof Necrosis) { // is current "cell" necrotic?
				necroticCell.add(obj);
			} 
		}
		
		if (necroticCell.size() > 0) {
			this.phagocytosis = phagocytosis - 1;
			for (int i = 0; i < necroticCell.size(); i++) {
				Object necrosisRandom = necroticCell.get(i);
				Necrosis.eatNecrosis(context, grid, necrosisRandom);
				this.numPhagocytosed++;
			}
		} else {
			m1Migrate();
		}
	}

	// M2 BEHAVIOR
	public void m2Behave(Context<Object> context) {
		// how do m2s migrate??
		// if m1 macs are gone but there is still necrosis - phagocytose
		//if (m2num > m1num && Necrosis.getNecrosis(context).size() > 0) {
			//this.m1Behave(context);
		//}
	}

	// MACROPHAGE PROLIFERATION
	public void macProliferate() {
		Context context = ContextUtils.getContext(this);
		GridPoint gridpt = grid.getLocation(this);
		Macrophage macrophage = new Macrophage(mcpSpatial, space, grid, this.phenotype, 5, this.phagocytosis,
				this.linked, this.buddy); // set phenotype = parent, age = 5; phagocytosed = same as parent; linked =
											// same as parent; buddy = same buddy as parent
		context.add(macrophage);
		// Object macrophage = this.clone(); // create copy of this agent (NetLogo
		// hatch)
		grid.moveTo(macrophage, gridpt.getX(), gridpt.getY());
	}

	// MACROPHAGE DEATH
	public void macDie() {
		Context context = ContextUtils.getContext(this);
		if (this.buddy != null) {
			// let the satellite buddy know this macrophage is gone - how?
		}
		context.remove(this);
	}

	// RETURN M1 AGENTS
	public static List<Object> getM1(Context<Object> context) { // Get a list of all the macrophages
		List<Object> m1Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 1) {
				m1Mac.add(obj);
			}
		}
		return m1Mac;
	}

	// RETURN M1 COUNT
	public int getM1Count() {
		return m1num;
	}

	// RETURN M2 AGENTS
	public static List<Object> getM2(Context<Object> context) { // Get a list of all the macrophages
		List<Object> m2Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 2) {
				m2Mac.add(obj);
			}
		}
		return m2Mac;
	}

	// RETURN M2 COUNT
	public int getM2Count() {
		return m2num;
	}

	// RETURN RESIDENT MAC AGENTS
	public static List<Object> getMres(Context<Object> context) { // Get a list of all the macrophages
		List<Object> mresMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 3) {
				mresMac.add(obj);
			}
		}
		return mresMac;
	}

	// RETURN RES MAC COUNT
	public int getMresCount() {
		return resMnum;
	}

	// RETURN M0 AGENTS (in the event that a phenotype was not assigned)
	public static List<Object> getM0(Context<Object> context) { // Get a list of all the macrophages
		List<Object> m0Mac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phenotype == 0) {
				m0Mac.add(obj);
			}
		}
		return m0Mac;
	}

	// RETURN MACROPHAGE PHENOTYPE
	public int getPhenotype() {
		return phenotype;
	}

	// RETURN APOP EATING MACROPHAGES
	public static List<Object> getM1ae(Context<Object> context) { // Get a list of all the macrophages
		List<Object> aeMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phagocytosis > 0) {
				aeMac.add(obj);
			}
		}
		return aeMac;
	}

	public double getM1aeCount() {
		Context context = ContextUtils.getContext(this);
		double m1ae = getM1ae(context).size();
		return m1ae;
	}

	// RETURN DEBRIS EATING MACROPHAGES
	public static List<Object> getM1de(Context<Object> context) { // Get a list of all the macrophages
		List<Object> deMac = new ArrayList<Object>(); // create a list of all the macrophage agents
		for (Object obj : context) {
			if (obj instanceof Macrophage && ((Macrophage) obj).phagocytosis < 0) {
				deMac.add(obj);
			}
		}
		return deMac;
	}
	
	public int getNumPhagocytosed() {
		return numPhagocytosed;
	}

	public double getM1deCount() {
		Context context = ContextUtils.getContext(this);
		double m1de = getM1de(context).size();
		return m1de ;
	}

}
