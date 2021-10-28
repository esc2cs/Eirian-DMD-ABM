/**
 * 
 */
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

/**
 * @author kbc7vn
 *
 */
public class Neutrophil {

	// Kyle's code: neutrophils-own [phagocytosis speed stay taxisx taxisy apoptosed
	// age]

	// Neutrophil parameters
	private static Grid<Object> grid;
	private static ContinuousSpace<Object> space;
	private static BufferedGridValueLayer mcpSpatial;
	public int phagocytosis; // neutrophil phagocytosis: starts at 0, 1+ converting to apoptotic neutrophil
								// after enough phagocytosis
	public int age; // neutrophil recruitment time- counter
	public int apoptosed; // neutrophil apoptosis counter
	public int numPhagocytosed; //temp to test fap phagocytosis

	/*
	 * // get growth factor numbers private static double[] growthFactors =
	 * GrowthFactors.getGrowthFactors(); // Recruiting growth factors private static
	 * double mcp = growthFactors[14]; private static double il1 = growthFactors[7];
	 * private static double il8 = growthFactors[8]; private static double ccl4 =
	 * growthFactors[12]; private static double cxcl1 = growthFactors[10]; private
	 * static double cxcl2 = growthFactors[9]; private static double gcsf =
	 * growthFactors[20]; // Deterring growth factors private static double lipoxins
	 * = growthFactors[22]; private static double mmp12 = growthFactors[19]; private
	 * static double resolvins = growthFactors[23]; private static double pge2 =
	 * growthFactors[27]; private static double lactoferrins = growthFactors[16];
	 * private static double il6 = growthFactors[13]; private static double il10 =
	 * growthFactors[21];
	 */

	public double numApopNeutrophils;
	public static double numNeutrophils;
	public double totalNeutrophilCount;

	// Recruitment-promoting: mcp1 (14), il1 (7), il8 (8), ccl4 (12), cxcl1 (10),
	// cxcl2 (9)
	// Recruitment-deterring: lipoxins (22), mmp12 (19), resolvins (23), pge2 (27),
	// lactoferrins (16), il6 (13), il10 (21)

	public Neutrophil(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, ContinuousSpace<Object> space, int age, int phagocytosis, int apoptosed) {
		this.grid = grid;
		this.space = space;
		this.age = age;
		this.phagocytosis = phagocytosis;
		this.apoptosed = apoptosed;
		this.mcpSpatial = mcpSpatial;
		this.numPhagocytosed = 0; //temp, using to test fap phagocytosis
	}
	
	// Called from Fiber (otherwise will not check for recruitment when Neutrophils leave system
	//@ScheduledMethod(start = 1, interval = 1, priority = 1, pick = 1)
	public static void neutrophilRecruitment(Context<Object> context) {// , Grid<Object> grid) { // recruitment/deter seesaw
		//Context context = ContextUtils.getContext(this);
		
		double[] growthFactors = GrowthFactors.getGrowthFactors();
		// Recruiting growth factors
		double mcp = growthFactors[14];
		double il1 = growthFactors[7];
		double il8 = growthFactors[8];
		double ccl4 = growthFactors[12];
		double cxcl1 = growthFactors[10];
		double cxcl2 = growthFactors[9];
		double gcsf = growthFactors[20];
		// Deterring growth factors
		double lipoxins = growthFactors[22];
		double mmp12 = growthFactors[19];
		double resolvins = growthFactors[23];
		double pge2 = growthFactors[27];
		double lactoferrins = growthFactors[16];
		double il6 = growthFactors[13];
		double il10 = growthFactors[21];

		// determine how strong the recruitment force is
		
		double recruit = il1 + gcsf + mcp + ccl4 + cxcl2 + cxcl1 + il8 + Necrosis.getPercentNecrotic(context)*100;
		double deter = lipoxins + resolvins + mmp12 + lactoferrins + pge2 + il10 + il6 / 2;
		double differential = recruit - deter/2;
		int addNeutrophils = 0;
		numNeutrophils = getNeutrophils(context).size();
		
		
		List<Object> necrosis = Necrosis.getNecrosis(context);
		List<Necrosis> newNecrosis = Necrosis.getNewNecrosis(context);
		if (differential > 0) {
			int max = Fiber.origFiberNumber * 2;
			int neutro_prob = RandomHelper.nextIntFromTo(0, (int) ((differential / 2) + (differential / 4) - 1)/25);
			if (neutro_prob > 10) {
				neutro_prob = (10) + RandomHelper.nextIntFromTo(0, 2);
			} else if (differential < 0) {
				if (Necrosis.getRecentPercentNecrotic(context) > 0 && numNeutrophils <= 0) {
				// if there was recent damage but neutrophils == 0 --> for instance at chronic damage-- then add
				addNeutrophils = (int) (Necrosis.getRecentPercentNecrotic(context) * 100 + Macrophage.getMres(context).size());
				neutro_prob = addNeutrophils;
				}
				System.out.println(neutro_prob);
			}
			
			//System.out.println("neutrophil: " + neutro_prob);
			while (neutro_prob > 0 && numNeutrophils < max) { // make some neutrophils
				// create_neutrophil();
				// find areas of necrosis to place neutrophil

				if (Necrosis.getNewNecrosis(context).size() > 0) {
					int index = RandomHelper.nextIntFromTo(0, newNecrosis.size() - 1);
					Object randomNecrosis = necrosis.get(index); // get the necrosis based on the random number chosen index
					
				
					GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location
					if (pt != null) {
						Neutrophil newNeutrophil = new Neutrophil(mcpSpatial, grid, space, 0, 0, 0); // change to neutrophil
						context.add((Neutrophil) newNeutrophil);
						grid.moveTo(newNeutrophil, pt.getX(), pt.getY());						
					} 
				} else if (necrosis.size() != 0) {
					int index = RandomHelper.nextIntFromTo(0, necrosis.size() - 1);
					Object randomNecrosis = necrosis.get(index); // get the necrosis based on the random number chosen index
					
				
					GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location
					if (pt != null) {
						Neutrophil newNeutrophil = new Neutrophil(mcpSpatial, grid, space, 0, 0, 0); // change to neutrophil
						context.add((Neutrophil) newNeutrophil);
						grid.moveTo(newNeutrophil, pt.getX(), pt.getY());						
					} 

				}
				neutro_prob = neutro_prob - 1;
				numNeutrophils = getNeutrophils(context).size();
			}
		}
	}

	@ScheduledMethod(start = 2, interval = 1, priority = 2)
	public void step() {
		Context context = ContextUtils.getContext(this);

		// Age neutrophil
		this.age = age + 1;

		if (this.apoptosed == 0) {
			// Find areas with necrosis

			// Find necrosis in this grid point
			MooreQuery<Object> query = new MooreQuery<Object>(grid, this, 1, 1); // get the location this neutrophil is
			// currently in
			Iterable<Object> iter = query.query(); // query the list of agents
			List<Object> necroticCell = new ArrayList<Object>();
			List<Object> fiberNgh = new ArrayList<Object>();
			List<Object> m1macNgh = new ArrayList<Object>();
			for (Object obj : iter) {
				if (obj instanceof Fiber) { // is neighbor a fiber?
					fiberNgh.add(obj);
				} else if (obj instanceof Macrophage) { // is neighbor a fiber?
					if (((Macrophage) obj).getPhenotype() == 1)
					m1macNgh.add(obj);
				}
			}
			
			GridPoint ptN = grid.getLocation(this);
			//get the objects on the same location as the neutrophil
			if (ptN != null) {
				Iterable<Object> iterPt = grid.getObjectsAt(ptN.getX(), ptN.getY());
				for (Object obj : iterPt) {
					if (obj instanceof Necrosis) { // is current "cell" necrotic?
						necroticCell.add(obj);
					} 
				}
			}


			if (necroticCell.size() > 0) {
				//int index = RandomHelper.nextIntFromTo(0, necroticCell.size() - 1);
				//Object necrosisRemove = necroticCell.get(index);
				//GridPoint pt = grid.getLocation(necrosisRemove);
				// neutrophil_phagocytose(context, pt);
				//Necrosis.eatNecrosis(context, grid, necrosisRemove); // should i add secondary necrosis instead of ECM?
				for (int i = 0; i < necroticCell.size(); i++) {
					Object necrosisRemove = necroticCell.get(i);
					Necrosis.eatNecrosis(context, grid, necrosisRemove);
					this.numPhagocytosed++;
					GridPoint pt = grid.getLocation(necrosisRemove);
					// secrete mcp when phagocytosing necrosis
					//mcpSpatial.set(5 + mcpSpatial.get(pt.getX(), pt.getY()), pt.getX(), pt.getY());		
				}
				if (fiberNgh.size() > 0 && m1macNgh.size() <= 0) {
					int index = RandomHelper.nextIntFromTo(0,  fiberNgh.size() - 1);
  					Fiber fiberRandom = (Fiber) fiberNgh.get(index);
					Necrosis newNecrosis = new Necrosis(mcpSpatial, grid, 1, 0); // Create a new Necrotic element, marked 'secondary'
					GridPoint pt = grid.getLocation(fiberRandom); // get location
					context.remove(fiberRandom); // remove the fiber and change to secondary necrosis
					context.add(newNecrosis); // Add new necrosis
					grid.moveTo(newNecrosis, pt.getX(), pt.getY()); 
				}
				this.phagocytosis = phagocytosis + 1;
				// secrete mcp when phagocytosing necrosis
				GridPoint pt = grid.getLocation(this);
				mcpSpatial.set(17 + mcpSpatial.get(pt.getX(), pt.getY()), pt.getX(), pt.getY());		
			} else {
				this.neutrophil_migrate();
			}

			int index = RandomHelper.nextIntFromTo(0, 8); // subtract 1 to account for Netlogo to Repast (3+3-1)
			if (index <= phagocytosis) { // After enough phagocytoses, it either becomes apoptotic or leaves
				this.apoptosed = apoptosed + 1;
				this.age = 0;
				// there is a 1/3 chance the neutrophil will now die
				int chance_die = RandomHelper.nextIntFromTo(0, 3);
				if (chance_die == 1) {
					context.remove(this);
				}
			}
		} else { // This is when apoptosed = 1, the promote M1 macrophages and deter neutrophils
			int chance_die = RandomHelper.nextIntFromTo(0, 3);
			if (chance_die == 1) {
				context.remove(this);
			}
			
		}

		// Once a neutrophil is 8 hours old, it has a chance of leaving, more likely the
		// longer its around
		if (this.age > 3 + 5 * 2) {
			int chance_die = RandomHelper.nextIntFromTo(0, 4);
			if (chance_die == 1) {
				context.remove(this);
			}
		}
		numNeutrophils = getNeutrophils(context).size();
	}

	// ADD NEW NEUTROPHILS TO NECROTIC TISSUE
	public void create_neutrophil() {
		Context context = ContextUtils.getContext(this);

		// find areas of necrosis to place neutrophil
		List<Object> necrosis = Necrosis.getNecrosis(context);
		if (necrosis.size() != 0) {
			int index = RandomHelper.nextIntFromTo(0, necrosis.size() - 1); // instead of random number, should they
																			// appear at edge necrosis first?
			Object randomNecrosis = necrosis.get(index); // get the necrosis based on the random number chosen index
															// within the list
			GridPoint pt = grid.getLocation(randomNecrosis); // Get the necrosis location

			Neutrophil newNeutrophil = new Neutrophil(mcpSpatial, grid, space, 0, 0, 0); // change to neutrophil
			context.add(newNeutrophil);
			grid.moveTo(newNeutrophil, pt.getX(), pt.getY());
		}

	}

	public void neutrophil_migrate() {
		Context  context = ContextUtils.getContext(this);
		// At each step neutrophils should sense surroundings and move toward
		// necrosis/high collagen
		MooreQuery<Object> query = new MooreQuery(grid, this, 3, 3); // get neighbors in a wider area
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		List<Object> necrosisNeighbor = new ArrayList<Object>();
		List<Object> openNeighbor = new ArrayList<Object>();
		List<Object> highCollNeighbor = new ArrayList<Object>();
		List<Object> ECMNeighbor = new ArrayList<Object>();
		for (Object neighbor : iter) {
			if (neighbor instanceof Necrosis) {
				necrosisNeighbor.add(neighbor); // if any neighbors are necrotic, add to the list
			} else if (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() > 1) { // neighbors with collagen
				highCollNeighbor.add(neighbor);
			} else if (neighbor instanceof ECM) { // neighbors with collagen
				ECMNeighbor.add(neighbor);
			} else {
				openNeighbor.add(neighbor);
			}
		}
		// If there is necrosis move there
		if (necrosisNeighbor.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, necrosisNeighbor.size() - 1);
			Object randomNeighbor = necrosisNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		// If there is no necrosis, move towards an area with collagen
		else if (highCollNeighbor.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, highCollNeighbor.size() - 1);
			Object randomNeighbor = highCollNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		// Otherwise just pick a random direction to move
		else if (ECMNeighbor.size() > 0) {
			int index = RandomHelper.nextIntFromTo(0, ECMNeighbor.size() - 1);
			Object randomNeighbor = ECMNeighbor.get(index);
			GridPoint pt = grid.getLocation(randomNeighbor);
			grid.moveTo(this, pt.getX(), pt.getY());
		}
		else {
		}
		

	}

	public void neutrophil_phagocytose(Context<Object> context, GridPoint pt) {
		this.phagocytosis = phagocytosis + 1;
		// Remove necrotic element
		context.remove(pt);
		// Replace with ECM element
		ECM newECM = new ECM(mcpSpatial, grid, 0.1, 0, 0);
		context.add(newECM);
		// Move ECM to where the necrosis was removed
		grid.moveTo(newECM, pt.getX(), pt.getY());
	}

	// Get the total apoptotic neutrophil count
	public static List<Object> getApoptosed(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		List<Object> apopNeutrophils = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof Neutrophil && ((Neutrophil) obj).apoptosed >= 1) {
				apopNeutrophils.add(obj);
			}
		}
		return apopNeutrophils;
	}
	
	public double getApoptosedCount() {
		Context context = ContextUtils.getContext(this);
		numApopNeutrophils = getApoptosed(context).size();
		return numApopNeutrophils;
	}

	// Get the total (non-apoptotic) neutrophil count
	public static List<Object> getNeutrophils(Context<Object> context) {
		//Context context = ContextUtils.getContext(this);
		List<Object> neutrophils = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof Neutrophil && ((Neutrophil) obj).apoptosed == 0) {
				neutrophils.add(obj);
			}
		}
		return neutrophils;
	}
	
	public double getNeutrophilCount() {
		return numNeutrophils;
	}
	
	/*
	public double getTotalNeutrophilCount() {
		totalNeutrophilCount = getNeutrophilCount() + getApoptosedCount();
		return totalNeutrophilCount;
	}*/

	public double getApoptosed() {
		return apoptosed;
	}

	public double getAge() {
		return age;
	}

	public double getPhagocytosis() {
		return phagocytosis;
	}
	
	public double getNumPhagocytosed() {
		return numPhagocytosed;
	}

}