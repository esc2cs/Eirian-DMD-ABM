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
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.query.space.grid.VNQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.SpatialMath;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;
import repast.simphony.valueLayer.GridValueLayer;

/**
 * @author Kelley Virgilio Satellite stem cells
 */
public class SSC {

	// SSC parameters
	private static Grid<Object> grid;
    private static BufferedGridValueLayer mcpSpatial; // MCP value layer	
    
    // activation signals/pressures:
	public static double sscActivation; // ssc activation pressure
	public static double sscMigration; // ssc migration pressure
	public static double sscProliferation; // ssc proliferation pressure
	public static double sscDifferentiation; // ssc differentiation pressure
	public static double sscQuiescence; // ssc quiescensce pressure
	public static double sscNiche; // tracks status of ssc niche
    public static double sscRecruitSatTemp; // saturation level for ssc recruitment

    // characteristics:
    private int active; // 0 = quiescent; 1-3 = activating SC; 4 = active + MyoD positive (myoblast)
    private int committed; // represents committed myogenic precursors 1= myoblasts in progress; 2 = myocyte in progress
	private int differentiated; // 0 = not differentiated; 1 = myogenin+; 2 = myocyte-- has differentiated for set amount of time; 3 = fused myocyte that can begin regenerating the fiber
    private int myf5neg; // marker for satellite cells that are Myf5-, do not differentiate, if equals 9 --> does not differentiate and restores the pool
	private int daughter; // tracks how many divisions 0- has not divided, 1 it has divided 1 time, 2 = 2 times
	private int migrated; // marks if the ssc has migrated or proliferated
	private int proteinAdd; // tracks how much protein (fiber elems) a ssc can add
    private int senescent; // tracks if the ssc is senescent
    public static int divideCount = 0; // count number of divisions

    // repair related:
    private int onEdge; // notes when a SSC has found a fiber edge in need of repair
	private int fiberNeedsRep; // attached to a fiber that needs repair
	private int fiberAssoc; // tracks what fiber the ssc is associated with so it can move when the fiber is regrowing

    // time tracking:
    private int migrationTime; // tracks how long the ssc has been migrating to fiber
	private int divideTime; // tracks how long it takes to divide
	private int diffTime; // tracks how long the ssc has been differentiated and fused
	private int sisterAssoc; // tracks how long sister cells should sit next to each other
    
    // time parameters (hours)
    private static final int timeToActive = 8; // hours it takes to become active [Cooper99 3-12 hours]
    private static final int timeToMigrate = 8; // hours it takes to migrate; with normal collagen density --> 6 + 2*1 = 8 (same as timeToActive) [Schuetz85 15 hours]
	private static final int timeToDivide = 10; // hours it takes to divide [Siegel11, Wang14]
	private static final int timeToDiff = 18; // hours it takes to differentiate [Flamini18 <24hrs]
    
    // other:
    private static final int maxProteinAdd = 60;//(int) (60 * (1 / Fiber.pax7Scale)); // This value should scale based on the size of the elements
	public static final double sscScale = 0.25;// .25 control

	public SSC(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, int active, int differentiated, int onEdge, int divideTime, int diffTime, int myf5neg,
			int daughter, int sisterAssoc, int proteinAdd, int fiberNeedsRep, int committed, int fiberAssoc,
			int senescent, int migrationTime) {
        this.grid = grid;
        this.mcpSpatial = mcpSpatial; 
		this.active = active; 
		this.differentiated = differentiated; 
		this.onEdge = onEdge; 
		this.divideTime = divideTime; 
		this.diffTime = diffTime; 
		this.myf5neg = myf5neg; 
		this.daughter = daughter; 
		this.sisterAssoc = sisterAssoc; 
		this.proteinAdd = proteinAdd; 
		this.fiberNeedsRep = fiberNeedsRep; 
		this.committed = committed; 
		this.fiberAssoc = fiberAssoc; 
		this.senescent = senescent; 
		this.migrationTime = migrationTime; 
	}

	// SSC BEHAVIORS at each time step/for each ssc agent
	@ScheduledMethod(start = 2, interval = 1)
	public void sscStep() {
		Context context = ContextUtils.getContext(this);
	    mcpSpatial = (BufferedGridValueLayer) context.getValueLayer("MCP Layer");

        // Active SSCs secrete mcp if the environment is inflammatory
        if (this.getActive() >= timeToActive && this.senescent == 0) {
            double inflamWeight = GrowthFactors.inflamWeight*100; //check if environment is inflammatory
            int randomInt = RandomHelper.nextIntFromTo(0, 100);
            if (randomInt <= inflamWeight) { // chance of secreting MCP increases with more inflammation
                GridPoint pt = grid.getLocation(this);
                mcpSpatial.set(10 + mcpSpatial.get(pt.getX(), pt.getY()), pt.getX(), pt.getY());
                double mcpHere = mcpSpatial.get(pt.getX(), pt.getY());
            }
        }

        // MOVE 
            // If active, not on edge, and not senescent
		if (this.getActive() >= timeToActive && this.getOnEdge() == 0 && this.senescent == 0) { // migration is tracked elsewhere
			move();
		}
        // MIGRATION 
            // Migration chance called from Fiber
		double timeToMigrateTemp = timeToMigrate - 2 + 2 * ECM.collagenDensity; // adjust ssc migration time based on the collagen density
		if (this.getMigrationTime() >= 1 && this.getMigrationTime() < timeToMigrateTemp) {
			this.setMigrationTime(this.getMigrationTime() + 1);
		}
		if (this.getMigrationTime() >= timeToMigrateTemp) {
			this.setActive(timeToActive); // set activation to 1
		}
        // SENSE ENVIRONMENT 
            // Determine if SSC is on a fiber that needs repair-- if so needs repair != -1
		int sscCountFiber = checkLocalEnviro();
		// SSC ACTIVATION
		if (this.getActive() == 0 && this.getMigrationTime() == 0) { // if it is migrating do not activate this way
			sscActivation(); // activates SSC if it is on fiber that needs repair or if pressure is high
		}
		if (this.getActive() >= 1 && this.getActive() < timeToActive) {
			setActive(this.getActive() + 1); // 1-3 = activating SC; 4 = active + MyoD positive (myoblast)
		}
        // SSC DIVISION 
            // Fully active, on edge of fiber
		if (this.getActive() >= timeToActive && this.senescent == 0) {
			sscDivision(sscCountFiber); 
		}
		if (this.getDivideTime() >= 1) {
			setDivideTime(this.getDivideTime() + 1); // to keep track of SSC division
		}
		// SISTER CELL ASSOCIATION
		if (this.getSisterAssoc() >= 1 && this.getSisterAssoc() < 8) {
			this.setSisterAssoc(getSisterAssoc() + 1); // count the hours until the sisters are not associated and can move/proliferate again
		}
		if (this.getSisterAssoc() >= 8) {
			this.setSisterAssoc(0); // reset to 0 after 8 hours to allow to move
		}
		// SSC DIFFERENTIATION
		    // Active; NOT dividing; NOT differentiated to myocyte
		if (this.getActive() >= timeToActive && this.getDivideTime() == 0 && this.getDiff() == 0 && this.senescent == 0) {
			sscDifferentiation();
		}
		if (this.getDiffTime() >= 1) {
			setDiffTime(this.getDiffTime() + 1); // to keep track of SSC differentiation
		}
		if (this.getDiffTime() > timeToDiff) { // once the set amount of differentiation time has occurred --> reset differentiation time to 0
			this.setDiffTime(0);
			// SSC --> MYOBLAST
			if (this.getCommitted() == 1) { // 1 is committing
				this.setCommitted(2); // 2 is a committed myogenic cell
			} else if (this.getCommitted() == 2) { // already a committed myogenic cell --> myocyte
				this.setDiff(2); // set differentiation to myocyte and make sure it is on an edge in order to fuse and regrow
			}
		}
		if (this.getOnEdge() == 1 && this.getDiff() == 2) { // if the ssc is on an edge and differentiated --> fused myocyte --> myofiber
			this.setDiff(3); // fused myocyte ready to regrow muscle
		}
		// FIBER REPAIR
		if (this.getDiff() == 3 && this.getProteinAdd() < maxProteinAdd && this.senescent == 0) { // if the ssc is differentiated- begin repairing the fiber
			Fiber.setFibersRepairing(1);
			// Disease state analysis:
			if (Fiber.macDepletion == 1 && ECM.collagenDensity > 1.5 && RandomHelper.nextIntFromTo(0, 1) < 1) { 
                // don't repair-- instead of turning off repair signal permanently
			} else {
				sscFiberRepair();
			}
		} else if (this.getProteinAdd() >= maxProteinAdd) {
			context.remove(this);
		}
		// QUIESCENCE
		    // Active, Not differentiated
		if (this.getActive() >= timeToActive && this.getDiff() == 0) {
			sscQuiescence();
		}
		// SSC and MYOBLAST APOPTOSIS
		    // Active, not differentiated ssc and myoblasts apoptose; not differentiating or dividing; and NOT on a fiber that still needs repair -- assume protected by fiber
		if (this.getActive() >= timeToActive && this.getDiff() == 0 && this.getDiffTime() == 0
				&& this.getDivideTime() == 0 && this.getFiberNeedsRep() != 1) {
			sscApoptosis();
		}
		// MIGRATE AWAY- RESTORE
		    // if SC have been here a long time and are Myf- then they have a chance of migrating back to restore populations on the rest of the cell
		if (this.getFiberNeedsRep() != 1 && this.getOnEdge() == 0 && this.getDiff() == 0 && this.getMyf5neg() == 9 && (RandomHelper.nextIntFromTo(0, 30) < 1)) {
			context.remove(this);
		}
	}

	// SSC SENSE ENVIRONMENT
	public int checkLocalEnviro() { // Each SSC determines its behavior in the context
		Context context = ContextUtils.getContext(this);
		double[] growthFactors = GrowthFactors.getGrowthFactors();
        // SSC PRESSURES
        double hgf = growthFactors[17];
        double vegf = growthFactors[18];
        double igf = growthFactors[2];
        double il6 = growthFactors[13];
        double tgfb = growthFactors[0];
        double pdgf = growthFactors[3];
        double gcsf = growthFactors[20];
        double il10 = growthFactors[21];
        double tnf = growthFactors[1];
        double fgf = growthFactors[29];
        double il4 = growthFactors[30];
        double il1 = growthFactors[7];
        double ecmProteins = growthFactors[6]; //fibronectin
        double mmp = growthFactors[4];
        double ifn = growthFactors[15];



		double activeTGF = GrowthFactors.getActiveTgf(InflamCell.getTick()); //is there a better way to call this now?
		// DISEASE STATE PARAMETERS:
		int sscCountFiber = 0;
		// SSC NICHE- function of fibronectin
		sscNiche = growthFactors[6];
		// SSC ACTIVATION- function of damage
		double sscActivation = 2*hgf*Fiber.getTotalFiberNumber(context) + fgf + igf; // 2*hgf + fgf + igf, hgf is scaled by fiber number
		
		//sscActivation = InflamCell.inflamCells[9]; // already scaled by fiber number
		// SSC DIVISION/PROLIFERATION:
		    // Division is induced by inflammatory factors/fibroblast factors: igf, fgf, tnf, ifn, il6, vegf, pdgf, gcsf, mmp
		    // Division is inhibited by anti-inflammatory factors: il10, tgf
		    // Division is weighted towards fgf and igf which are required for cell cycle entry
		//sscProliferation = (3 * igf + 3 * fgf + tnf + ifn + il6 + vegf + pdgf + gcsf - il10 - 4 * activeTGF) / Fiber.origFiberNumber;
		sscProliferation = (3 * igf + 3 * fgf + tnf + ifn + il6 + vegf + pdgf + gcsf - il10 - 2 * activeTGF) / Fiber.origFiberNumber;
		// SSC CELL CYCLE EXIT/DIFFERENTIATION:
		    // Cell cycle exit is induced by il10, il4
		    // Cell cycle exit is blocked by fgf, igf, hgf, ifn, tnf
		        // Note: "blocked" factors do not necessarily imply that they block differentiation-- promotes ssc to stay in cell cycle which will still produce differentiated cells eg IGF promotes prolif and diff in lit, negative because it promotes cells to stay in cycle
		sscDifferentiation = (-2 * fgf - 2 * igf - 2 * hgf - ifn - tnf + 4 * il10 + il4) / Fiber.origFiberNumber;
		// SSC MIGRATION
		    // Migration is induced by hgf, igf, fgf, mmps
		    // Migration is blocked by tgf
		sscMigration = ((hgf + igf + fgf) - 2 * (activeTGF) + mmp / ECM.collagenDensity) / Fiber.origFiberNumber * .85;
		double fiberNorm = Fiber.origFiberNumber;
		if (Fiber.origFiberNumber < getActiveSSCs(context).size()) {
			// If there are more SSCs than fibers-- then normalize by the number of SSCs - check Reimann et al. paper
			fiberNorm = getActiveSSCs(context).size(); 
		}
		sscRecruitSatTemp = sscActivation * Fiber.origFiberNumber; // (fiberNorm * 3); 
		// SSC QUIESCENCE:
		    // Function of the amount of the muscle that is recovered
		sscQuiescence = (Fiber.origCsaMean * Fiber.origFiberNumber - Fiber.getFiberElems(context).size()) / (Fiber.origFiberNumber) * 50;
		if (sscQuiescence < 50) {
			sscQuiescence = 50;
		}
		// Check if the SSC is on a fiber that needs repair
		// Get the SSC neighbors and see what fiber they are attached to:
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get neighbors
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		int needsRepairTemp = -1;
		int sscCountOnFiberTemp = 0;
		int onEdgeTemp = -1;
		Object fiberAtSSC = null;
		for (Object obj : iter) { // go through the neighbors of the satellite stem cells and determine what fiber it is on
			if (obj instanceof Fiber && ((Fiber) obj).needsRepair == 1 && needsRepairTemp == -1) {
				this.setFiberNeedsRep(1);
				sscCountFiber = ((Fiber) obj).getsscCountOnFiber();
				needsRepairTemp = 1; // update value so it changes
			}
			if (obj instanceof Fiber && onEdgeTemp == -1) {
				onEdgeTemp = 1; // update
			}
		}
		// if none of the neighbors need repair
		if (needsRepairTemp == -1) {
			this.setFiberNeedsRep(0);
		}
		if (onEdgeTemp == -1) {
			// the ssc is not near ANY fiber
			if (this.fiberAssoc != 0) {
				// If it is associated with a fiber, move it there-- because of restructuring it might have been separated
				List<Object> fiberElems = Fiber.getElemInFiber(this.fiberAssoc, context);
				if (fiberElems.size() > 0) {
					int rando = RandomHelper.nextIntFromTo(0, fiberElems.size() - 1);
					Object randoFib = fiberElems.get(rando);
					if (randoFib != null) {
						grid.moveTo(this, grid.getLocation(randoFib).getX(), grid.getLocation(randoFib).getY()); // just pick a random spot and move there
					}
				}
			} else {
				setOnEdge(0);
			}
		}
		return sscCountFiber; // returns the number of SSCs on the fiber this SSC is nearest
	}

	// ACTIVATION
	public void sscActivation() {
		if (this.getFiberNeedsRep() == 1) { // SSC activation-- if the fiber needs repair-- set active to 1
			// SSC will activate if it is on a fiber that needs repair
			this.setActive(1); 
		}
		// Otherwise there is just a chance of activation because of nearby damage
		else if (sscActivation > Fiber.origFiberNumber * 3	&& RandomHelper.nextIntFromTo(0, (50 - (int) Math.floor(sscActivation))) < 1) { // SSC activation
			this.setActive(1); // activate SSC if pressure is high
		}
	}

	// SSC DIVISION
	public void sscDivision(int sscCountFiber) {
		Context context = ContextUtils.getContext(this);
        double[] growthFactors = GrowthFactors.getGrowthFactors();
        double fgf = growthFactors[29];
        double igf = growthFactors[2];
		
        //if (sscCountFiber < 6 && this.getOnEdge() == 1 && this.fiberNeedsRep == 1 && fgf > Fiber.origFiberNumber * 4 && igf > Fiber.origFiberNumber * 10) { // ADDED REQUIREMENT FOR IGF AND FGF
		if (sscCountFiber < 15 && this.getOnEdge() == 1 && this.fiberNeedsRep == 1 && fgf > Fiber.origFiberNumber * 4 && igf > Fiber.origFiberNumber * 10) { // ADDED REQUIREMENT FOR IGF AND FGF
			// Only requirement for this function is that the ssc is active--- in order to make a new division check that it is on an edge that needs repair
			if (this.getDivideTime() == 0) { // if active and not already proliferating
				// It takes 8-14 hours to activate originally - first division takes 18-24 hours while each subsequent division takes 10 hours
				int divNumberEffect = 1;
				if (this.daughter >= 1) {
					// if the daughter has already divided twice- change the chance of division
					divNumberEffect = (this.daughter + 1);
				}
				int divChance = (int) ((220 - sscProliferation)); // Chance of division weighted by division signal; 8302016 parameter analysis solution
				//if (divChance < 60) { // threshold for division chance- 8302016 parameter analysis solution
				//	divChance = 60;
				//}
				if (divChance < 80) { // threshold for division chance- trying new numbers now
					divChance = 80;
				}
				if (RandomHelper.nextIntFromTo(0,(int) (divChance * Math.sqrt(divNumberEffect) * Fiber.regenCapacity)) < 1 && this.getDiff() == 0) { // Chance of dividing
					this.setDivideTime(1); // set divide time
					this.setDaughter(1); // Should ssc that divided limit the number of divisions it can do
					divideCount++; // count number of divisions
				}
			}
			// Even if the fiber already is repaired, it can still complete its division
			else if (this.getDivideTime() >= timeToDivide) { // It takes on average 10 hours for the cell to divide
				this.setDivideTime(0); // reset proliferating to 0 so it can proliferate again
				// Get ECM location- proliferating SSC should be added to neighbor of SSC and ECM
				List<Object> neighbor = new ArrayList<Object>();
				VNQuery<Object> query2 = new VNQuery(grid, this, 1, 1); // get neighbors
				Iterable<Object> iter2 = query2.query(); // query the list of agents that are the neighbors
				// Choose a random neighbor and place the agent there
				for (Object neighborIter : iter2) {
					neighbor.add(neighborIter);
				}
				int index = RandomHelper.nextIntFromTo(0, neighbor.size() - 1); // draw random number with randomHelper
				Object randomNeigh = neighbor.get(index); // get the ecm based on the random number chosen-- index
															// within the list
				GridPoint pt = grid.getLocation(randomNeigh); // Get the ecm location
				// Fibronectin dependence- ssc symmetric division depends on the presence of
				// fibronectin
				double fibronectinFactor = 1.;
				if (sscNiche / Fiber.origFiberNumber > 12) {
					fibronectinFactor = 1.;
				} else {
					fibronectinFactor = (sscNiche / Fiber.origFiberNumber) / 20.;
					//System.out.println(fibronectinFactor);
				}
				// SC DAUGHTER FATE DETERMINATION
				// IF IT IS A SSC
				int chanceSymmetric = 90; // 10% chance of asymmetric division in DMD [Dumont2015 ]
				if (this.getCommitted() == 0) { // If it is a SSC (not committed)
					if (RandomHelper.nextIntFromTo(0, 100) > (chanceSymmetric) * fibronectinFactor && this.getMyf5neg() != 9) { // 10% of SC go through asymmetric division 90% go through symmetric division with a SC and committed daughter
						// MYOBLAST- ASYMMETRIC DIVISION
						SSC sscNewAsymm = new SSC(mcpSpatial, grid, 1, 0, 1, 0, 0, RandomHelper.nextIntFromTo(0, 8), (this.getDaughter() + 1), 1, 0, 0, 2, this.getFiberAssoc(), 0, 0); // adds an active SSC
						// cell is not differentiated, but committed, can't be myf5- (myf5- == 9)
						context.add((SSC) sscNewAsymm); // add to context
						grid.moveTo(sscNewAsymm, pt.getX(), pt.getY()); // add the new ssc to that location
						this.setSisterAssoc(1);
					} else {
						// SSC SYMMETRIC DIVISION
						SSC sscNewSymm = new SSC(mcpSpatial, grid, 1, 0, 1, 0, 0, RandomHelper.nextIntFromTo(0, 9), (this.getDaughter() + 1), 5, 0, 0, 0, this.getFiberAssoc(), 0, 0); // adds an active SSC
						context.add((SSC) sscNewSymm); // add to context
						grid.moveTo(sscNewSymm, pt.getX(), pt.getY()); // add the new ssc to that location
						this.setDaughter(this.getDaughter() + 1); // note that the current cell has divided
						sscNewSymm.setSisterAssoc(5); // only associates with a sister sc for 3 hours, so it starts at hour 5
						this.setSisterAssoc(5);
					}
				}
				// IF IT IS A MYOBLAST
				else if (this.getCommitted() == 2) {
					if (RandomHelper.nextIntFromTo(0, 100) > (chanceSymmetric) * fibronectinFactor) { // 10% of SSC go through asymmetric division. 90% go through symmetric division with a SSC and committed daughter
																												
						// MYOCYTE- ASYMMETRIC DIVISION
							// can't be myf5- (myf5- == 9)
						SSC sscNewAsymm = new SSC(mcpSpatial, grid, 1, 1, 1, 0, 1, RandomHelper.nextIntFromTo(0, 8), (this.getDaughter() + 1), 1, 0, 0, 2, this.getFiberAssoc(), 0, 0); // adds an active SSC
						context.add((SSC) sscNewAsymm); // add to context
						grid.moveTo(sscNewAsymm, pt.getX(), pt.getY()); // add the new ssc to that location
						this.setSisterAssoc(1);				

					} else {
						// MYOBLAST
						    // can't be myf5- (myf5- == 9)
						SSC sscNewSymm = new SSC(mcpSpatial, grid, 1, 0, 1, 0, 0, RandomHelper.nextIntFromTo(0, 8), (this.getDaughter() + 1), 5, 0, 0, 2, this.getFiberAssoc(), 0, 0); // adds an active SSC
						context.add((SSC) sscNewSymm); // add to context
						grid.moveTo(sscNewSymm, pt.getX(), pt.getY()); // add the new ssc to that location
						this.setDaughter(this.getDaughter() + 1); // note that the current cell has divided
						sscNewSymm.setSisterAssoc(5); // only associates with a sister sc for 3 hours, so it starts at hour 5
						this.setSisterAssoc(5);
					}
				}
			}
		}
	}

	// QUIESCENCE- If the SSC is on a fiber that does not need to be repaired there is a chance for it to become quiescent again
	public void sscQuiescence() {
		// if the SSC/myoblast is on a fiber that does not need to be repaired --> chance to return to quiescense
		Context context = ContextUtils.getContext(this);
		if (this.fiberNeedsRep != 1 && this.onEdge == 1 && RandomHelper.nextIntFromTo(0, 100) < 1) {
			this.setActive(0); // if there is a strong signal from repaired fibers than the SSC has a chance to become quiescent
		}
		// if the SC is not on a fiber edge but there is a large quiescence signal --> chance of quiescence
		else if (this.onEdge == 0 && this.fiberNeedsRep == 0
				&& RandomHelper.nextIntFromTo(0, (int) (sscQuiescence)) < 1) {
			// include fiber does not need repair in case it is one off edge but still associated with a fiber repair
			this.setActive(0);
		}
	}

	// APOPTOSIS
	public void sscApoptosis() {
		// ssc and myoblasts have a chance of apoptosing- M1s stop apoptosis if there are more myoblasts than M1 cells --> chance M1 macrophages protect SSC from apoptosis
		Context context = ContextUtils.getContext(this);
		double diffM1MyobTemp = (getMyoblastCount() - (InflamCell.inflamCells[3] + InflamCell.inflamCells[4] + InflamCell.inflamCells[5] + GrowthFactors.m1MacAdded)) / (getMyoblastCount() + 1);
		int apopChanceTemp = (int) (250 * (1 - diffM1MyobTemp));
		if (diffM1MyobTemp > 0 && RandomHelper.nextIntFromTo(0, apopChanceTemp) < 1) {
			context.remove(this);
		}
	}

	// EXIT CELL CYCLE, TERMINALLY DIFFERENTIATE
	public void sscDifferentiation() {
		// Two types of differentiation: SSC --> myoblast and myoblast --> myocyte
		int diffChanceTemp = (int) (0 - sscDifferentiation); // chance of differentiation weighted by differentiation signal
		if (diffChanceTemp < 5) { // low limit on chance of differentiation
			diffChanceTemp = 5;
		}
		if (RandomHelper.nextIntFromTo(0, (int) (diffChanceTemp * (Fiber.pax7Scale))) < 1 && getMyf5neg() != 9) {
			// 10% of population will not differentiate (can still result in a Myf5+ daughter cell during division though
			// SSC --> MYOBLAST
			if (this.getCommitted() == 0) {
				this.setCommitted(1);
				this.setDiffTime(1);
			}
			// MYOBLAST --> MYOCYTE
			else if (this.getCommitted() == 2) {
				this.setDiff(1); // will become a myocyte once it has fully differentiated
				this.setDiffTime(1);
			}
		}
	}

	// Called from Fiber class with 'pick 1'
	public static void sscMigrationChance(Context<Object> context) { // SSC MIGRATION
		int migrChanceTemp = (int) (150 - 2 * sscMigration); // chance of ssc migration based on migration signal // 8302016 parameter analysis result
		if (migrChanceTemp < 30) { // lower threshold for chance - 8302016 parameter analysis result
			migrChanceTemp = 30;
		}
		// saturation count = active SSCs + ssc migrating for 8 hours (do not want slowed migration to increase number recruited count number of recruited sscs
		int sscCountMigrated = 0;
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).active == 0 && ((SSC) obj).migrationTime >= timeToMigrate) { // number of SSCs that are migrating but delayed and have not arrived
				sscCountMigrated = sscCountMigrated + 1;
			}
		}
		int sscCountFiber = getActiveSSCs(context).size() + sscCountMigrated;
		if (sscMigration > 0 && RandomHelper.nextIntFromTo(0, migrChanceTemp) < 1 && sscCountFiber < sscRecruitSatTemp * Fiber.pax7Scale) { // decrease sscActivation threshold for low damage
			// Get a list of all the ECM that borders a fiber that needs repair
			List<Object> fibers = Fiber.getFiberElems(context); // Get a list of all the fibers
			List<Object> fiberBorderDamaged = new ArrayList<Object>();
			List<Object> fiberBorderDamaged1 = new ArrayList<Object>();
			List<Object> fiberBorderDamaged2 = new ArrayList<Object>();
			List<Object> fiberBorderDamaged3 = new ArrayList<Object>();
			for (Object fiberDamCheck : fibers) {
				// Preferentially recruited to fibers with low ssc counts and high damage
				if (((Fiber) fiberDamCheck).getBorder() == 1 && ((Fiber) fiberDamCheck).getNeedsRepair() == 1 && ((Fiber) fiberDamCheck).getsscCountOnFiber() < 1) {
					// if it is a fiber edge and it needs repair- add to the list
					fiberBorderDamaged.add(fiberDamCheck);
				}
				else if (((Fiber) fiberDamCheck).getBorder() == 1 && ((Fiber) fiberDamCheck).getNeedsRepair() == 1 && ((Fiber) fiberDamCheck).getsscCountOnFiber() == 1) {
					// if it is a fiber edge and it needs repair- add to the list
					fiberBorderDamaged1.add(fiberDamCheck);
				}
				else if (((Fiber) fiberDamCheck).getBorder() == 1 && ((Fiber) fiberDamCheck).getNeedsRepair() == 1 && ((Fiber) fiberDamCheck).getsscCountOnFiber() <= 3) {
					// if it is a fiber edge and it needs repair- add to the list
					fiberBorderDamaged2.add(fiberDamCheck);
				}
				else if (((Fiber) fiberDamCheck).getBorder() == 1 && ((Fiber) fiberDamCheck).getNeedsRepair() == 1 && ((Fiber) fiberDamCheck).getsscCountOnFiber() < 6) {
					// Set limit on number of SSC that would migrate to a fiber if it is a fiber edge and it needs repair- add to the list
					fiberBorderDamaged3.add(fiberDamCheck);
				}
			}
			// Determination of where to migrate to:
			    // Migrates along the fiber to an edge with damage
			    // Choose a random damaged fiber edge and add a SC in the nearby ECM
			    // Only adds a SC if there is a need; Only migrates if there is not a saturatednumber of sscs
			if (fiberBorderDamaged.size() > 0) {
				SSC sscNew = new SSC(mcpSpatial, grid, 0, 0, 0, 0, 0, RandomHelper.nextIntFromTo(0, 9), 0, 0, 0, 0, 0, 0, 0, 1); // adds a quiescent SSC that has not migrated
				// Chance is dependent on the number of ssc
				context.add(sscNew); // add to context
				int randomInt = RandomHelper.nextIntFromTo(0, fiberBorderDamaged.size() - 1);
				Object damBorderRand = fiberBorderDamaged.get(randomInt); // get the random damaged border fiber
				((SSC) sscNew).setFiberAssoc(((Fiber) damBorderRand).getFiberNumber());
				// Get the ecm neighbor to place the SSC there
				MooreQuery<Object> query = new MooreQuery(grid, damBorderRand, 1, 1); // get neighbors
				Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
				List<Object> openNeighbors = new ArrayList<Object>();
				for (Object neighbors : iter) {
					if (neighbors instanceof ECM || neighbors instanceof Necrosis) {
						openNeighbors.add(neighbors);
					}
				}
				int int2 = RandomHelper.nextIntFromTo(0, openNeighbors.size() - 1);
				if (openNeighbors.size() > 0) {
					Object ecmNearDam = openNeighbors.get(int2); // get the ecm near the damaged fiber
					GridPoint ptECM = grid.getLocation(ecmNearDam); // Get the ecm location
					grid.moveTo(sscNew, ptECM.getX(), ptECM.getY()); // add the new ssc to that location
					sscNew.setMigrated(1);
				}
			} else if (fiberBorderDamaged1.size() > 0) {
				SSC sscNew = new SSC(mcpSpatial, grid, 0, 0, 0, 0, 0, RandomHelper.nextIntFromTo(0, 9), 0, 0, 0, 0, 0, 0, 0, 1); // adds a quiescent SSC that has not migrated
				// Chance is dependent on the number of ssc
				context.add(sscNew); // add to context
				int randomInt = RandomHelper.nextIntFromTo(0, fiberBorderDamaged1.size() - 1);
				Object damBorderRand1 = fiberBorderDamaged1.get(randomInt); // get the random damaged border fiber
				((SSC) sscNew).setFiberAssoc(((Fiber) damBorderRand1).getFiberNumber());
				// Get the ecm neighbor to place the SSC there
				MooreQuery<Object> query1 = new MooreQuery(grid, damBorderRand1, 1, 1); // get neighbors
				Iterable<Object> iter1 = query1.query(); // query the list of agents that are the neighbors
				List<Object> openNeighbors1 = new ArrayList<Object>();
				for (Object neighbors1 : iter1) {
					if (neighbors1 instanceof ECM || neighbors1 instanceof Necrosis) {
						openNeighbors1.add(neighbors1);
					}
				}
				int int1 = RandomHelper.nextIntFromTo(0, openNeighbors1.size() - 1);
				Object ecmNearDam = openNeighbors1.get(int1); // get the ecm near the damaged fiber
				GridPoint ptECM1 = grid.getLocation(ecmNearDam); // Get the ecm location
				grid.moveTo(sscNew, ptECM1.getX(), ptECM1.getY()); // add the new ssc to that location
				sscNew.setMigrated(1);
			} else if (fiberBorderDamaged2.size() > 0 && RandomHelper.nextIntFromTo(0, 5) < 1) { // if all the damaged fibers have SSC on them than pick a random one
				// Chance is dependent on the number of ssc
				SSC sscNew = new SSC(mcpSpatial, grid, 0, 0, 0, 0, 0, RandomHelper.nextIntFromTo(0, 9), 0, 0, 0, 0, 0, 0, 0, 1); // adds a quiescent SSC that has not migrated 
				context.add(sscNew); // add to context
				int randomInt = RandomHelper.nextIntFromTo(0, fiberBorderDamaged2.size() - 1);
				Object damBorderRand2 = fiberBorderDamaged2.get(randomInt); // get the random damaged border fiber
				((SSC) sscNew).setFiberAssoc(((Fiber) damBorderRand2).getFiberNumber());
				// Get the ecm neighbor to place the SSC there
				MooreQuery<Object> query2 = new MooreQuery(grid, damBorderRand2, 1, 1); // get neighbors
				Iterable<Object> iter2 = query2.query(); // query the list of agents that are the neighbors
				List<Object> openNeighbors2 = new ArrayList<Object>();
				for (Object neighbors2 : iter2) {
					if (neighbors2 instanceof ECM || neighbors2 instanceof Necrosis) {
						openNeighbors2.add(neighbors2);
					}
				}
				int int2 = RandomHelper.nextIntFromTo(0, openNeighbors2.size() - 1);
				Object ecmNearDam2 = openNeighbors2.get(int2); // get the ecm near the damaged fiber
				GridPoint ptECM2 = grid.getLocation(ecmNearDam2); // Get the ecm location
				grid.moveTo(sscNew, ptECM2.getX(), ptECM2.getY()); // add the new ssc to that location
				sscNew.setMigrated(1);
			} else if (fiberBorderDamaged3.size() > 0 && RandomHelper.nextIntFromTo(0, 10) < 1) { // if all the damaged fibers have SSC on them, then pick a random one
				// Chance is dependent on the number of ssc
				SSC sscNew = new SSC(mcpSpatial, grid, 0, 0, 0, 0, 0, RandomHelper.nextIntFromTo(0, 9), 0, 0, 0, 0, 0, 0, 0, 1); // adds a quiescent SSC that has not migrated
				context.add(sscNew); // add to context
				int randomInt = RandomHelper.nextIntFromTo(0, fiberBorderDamaged3.size() - 1);
				Object damBorderRand3 = fiberBorderDamaged3.get(randomInt); // get the random damaged border fiber
				((SSC) sscNew).setFiberAssoc(((Fiber) damBorderRand3).getFiberNumber());
				// Get the ecm neighbor to place th SSC there
				MooreQuery<Object> query3 = new MooreQuery(grid, damBorderRand3, 1, 1); // get neighbors
				Iterable<Object> iter3 = query3.query(); // query the list of agents that are the neighbors
				List<Object> openNeighbors3 = new ArrayList<Object>();
				for (Object neighbors3 : iter3) {
					if (neighbors3 instanceof ECM || neighbors3 instanceof Necrosis) {
						openNeighbors3.add(neighbors3);
					}
				} 
				if (openNeighbors3.size() > 0) {
					int int2 = RandomHelper.nextIntFromTo(0, openNeighbors3.size() - 1);
					Object ecmNearDam3 = openNeighbors3.get(int2); // get the ecm near the damaged fiber
					GridPoint ptECM3 = grid.getLocation(ecmNearDam3); // Get the ecm location
					grid.moveTo(sscNew, ptECM3.getX(), ptECM3.getY()); // add the new ssc to that location
					sscNew.setMigrated(1);
				}
			}
		}
	}

	// Differentiated myocytes repair muscle
	public void sscFiberRepair() {
		// Determine what fiber it is on and regrow
		if (RandomHelper.nextIntFromTo(0, (int) (4 * Fiber.pax7Scale)) < 1) {
			Context context = ContextUtils.getContext(this);
			MooreQuery<Object> query = new MooreQuery(grid, this, 3, 3); // get neighbors
			// Search a larger area because SSC division sometimes ends up wit a fiber slightly further away
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			List<Object> nearbyFiber = new ArrayList<Object>();
			int fiberNumber = -1;
			Object fiberAtSSC = null;
			for (Object obj : iter) { // go through the neighbors of the satellite stem cells
				if (obj instanceof Fiber && ((Fiber) obj).needsRepair == 1) {
					nearbyFiber.add((Fiber) obj);
				}
			}
			if (nearbyFiber.size() > 0) {
				// if there is an adjacent fiber add at that location
				int randomInt = RandomHelper.nextIntFromTo(0, nearbyFiber.size() - 1);
				fiberAtSSC = nearbyFiber.get(randomInt);
				fiberNumber = ((Fiber) fiberAtSSC).getFiberNumber();
				setFiberAssoc(fiberNumber);
				((Fiber) fiberAtSSC).addFiberElem(fiberNumber, context, 0); // add a fiber element/agent at this fiber
				this.setProteinAdd(getProteinAdd() + 1); // tracks how many fiber elements are added
			}
		}
	}

	public void move() { // SSC move towards fibers that need repair
		// SSC: NON-DIFFERENTIATED
		if (RandomHelper.nextIntFromTo(0, (int) (ECM.collagenDensity)) <= 1) { // Move IF collagen density is low
			Context context = ContextUtils.getContext(this);
			MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // get neighbors
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			MooreQuery<Object> queryW = new MooreQuery(grid, this, 10, 10); // get neighbors within a wider range
			Iterable<Object> iterW = queryW.query(); // query the list of agents that are the neighbors
			List<Object> fiberNeighborW = new ArrayList<Object>();
			List<Object> damageNeighborW = new ArrayList<Object>();
			List<Object> damageNeighbor = new ArrayList<Object>();
			List<Object> openNeighbor = new ArrayList<Object>();
			// Is it near a fiber that needs repair?
			double lowCollagen = 1.5;
			for (Object neighbor : iter) {
				//if (neighbor instanceof Fiber && ((Fiber) neighbor).getNeedsRepair() == 1 && ((Fiber) neighbor).getsscCountOnFiber() < 6) { // if it borders a fiber that needs repair
				if (neighbor instanceof Fiber && ((Fiber) neighbor).getNeedsRepair() == 1 && ((Fiber) neighbor).getsscCountOnFiber() < 15) { // if it borders a fiber that needs repair
					setOnEdge(1); // marks that the ssc is on a fiber edge
					List<Object> elemInFiberTemp = ((Fiber) neighbor).getElemInFiber(((Fiber) neighbor).getFiberNumber(), context); // get the elements in this fiber
					for (Object fiberElems : elemInFiberTemp) {
						// Only count if the SSC is NOT myf9 --> since this will not contribute to differentiation
						if (this.getMyf5neg() != 9) {
							// marks that the ssc is on the fiber edge--- can only do this one time
							((Fiber) fiberElems).setsscCountOnFiber(((Fiber) neighbor).getsscCountOnFiber() + 1); 
						}
					}
					return;
				} else {
					
					if (neighbor instanceof ECM) { // also get a list of ecm
						openNeighbor.add(neighbor); // Default: move within ecm
					}
					if (neighbor instanceof Necrosis || (neighbor instanceof ECM && ((ECM) neighbor).getCollagen() < lowCollagen)) { 
						// also get a list of necrotic and low collagen ECM neighbors in case
						damageNeighbor.add(neighbor); // Default: move within ecm
					}
				}
			}

			if (this.getOnEdge() == 0) { // if it is not on the fiber edge
				for (Object neighborW : iterW) {
					//if (neighborW instanceof Fiber && ((Fiber) neighborW).getNeedsRepair() == 1 && ((Fiber) neighborW).getsscCountOnFiber() < 6) {
					if (neighborW instanceof Fiber && ((Fiber) neighborW).getNeedsRepair() == 1 && ((Fiber) neighborW).getsscCountOnFiber() < 15) {
						fiberNeighborW.add(neighborW);
					}
					if (neighborW instanceof Necrosis || (neighborW instanceof ECM && ((ECM) neighborW).getCollagen() < lowCollagen)) {
						damageNeighborW.add(neighborW);
					}
				}
				if (fiberNeighborW.size() > 0) { // If there is a fiber within the wider range
					int indexW = RandomHelper.nextIntFromTo(0, fiberNeighborW.size() - 1);
					Object randomNeighborW = fiberNeighborW.get(indexW); // choose one of the fibers within the range
					GridPoint ptW = grid.getLocation(randomNeighborW); // get current location
					moveTowards(ptW); // Move ssc towards the fiber neighbor
				} else { // if there is not a fiber edge needing repair in the area then move towards necrosis
					// If there is no fiber needed repair close by-- chance of dying/leaving
					if (damageNeighbor.size() > 0) { // if there is necrosis move there
						int index = RandomHelper.nextIntFromTo(0, damageNeighbor.size() - 1);
						Object randomNeighbor = damageNeighbor.get(index);
						GridPoint pt = grid.getLocation(randomNeighbor);
						grid.moveTo(this, pt.getX(), pt.getY());
					} else if (damageNeighborW.size() > 0) { // If there is damage in the wider range
						int indexW = RandomHelper.nextIntFromTo(0, damageNeighborW.size() - 1);
						Object randomNeighborW = damageNeighborW.get(indexW); // choose one of the fibers within the range
						GridPoint ptW = grid.getLocation(randomNeighborW); // get current location
						moveTowards(ptW); // Move ssc towards the fiber neighbor
					} else if (openNeighbor.size() > 0) { // Otherwise just pick a random direction of ecm go to it
						int index = RandomHelper.nextIntFromTo(0, openNeighbor.size() - 1);
						Object randomNeighbor = openNeighbor.get(index);
						GridPoint pt = grid.getLocation(randomNeighbor);
						grid.moveTo(this, pt.getX(), pt.getY());
					}
				}
			}
		}
	}

	public void moveTowards(GridPoint pt) {
		// Need to make sure ssc move in ecm and necrosis only
		GridPoint myPoint = grid.getLocation(this);
		if ((pt.getY() - myPoint.getY()) >= 1) {
			try{
				grid.moveTo(this, myPoint.getX(), myPoint.getY() + 1);
			} catch (Exception e) {}
		} else if ((pt.getY() - myPoint.getY()) < 1) {
			try {
				grid.moveTo(this, myPoint.getX(), myPoint.getY() - 1);
			} catch (Exception e) {}
		}
		GridPoint myPointNew = grid.getLocation(this); // update location
		if ((pt.getX() - myPointNew.getX()) >= 1) {
			try {
				grid.moveTo(this, myPointNew.getX() + 1, myPointNew.getY());
			} catch (Exception e) {}
		} else if ((pt.getX() - myPointNew.getX()) < 1) {
			try {
				grid.moveTo(this, myPointNew.getX() - 1, myPointNew.getY());
			} catch (Exception e) {}
		}
	}

	public static void initialize(int origFiberNumber, Context<Object> context, Grid<Object> grid) {
		// Initialize ssc with non-active/quiescent cells that are attached to the muscle fiber
		for (int i = 0; i < Math.ceil(origFiberNumber * sscScale * Fiber.pax7Scale); i++) { // less ssc at quiescence than fibroblasts
			context.add(new SSC(mcpSpatial, grid, 0, 0, 0, 0, 0, RandomHelper.nextIntFromTo(0, 9), 0, 0, 0, 0, 0, 0, 0, 0)); // Add the set number of non-active ssc to context																									
		}
		List<Object> ecms = ECM.getECM(context); // Get a list of all the ecm agents to place cells on the ECM
		List<Object> ecmAtBorder = new ArrayList<Object>();
		List<Object> fiberNumberTemp = new ArrayList<Object>();
		for (Object ecmIter : ecms) {
			VNQuery<Object> query = new VNQuery(grid, ecmIter, 1, 1); // Find the 4 neighbors and determine if the ecm borders a fiber
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			for (Object neighborCheck : iter) { // go through the list of neighbors
				if (neighborCheck instanceof Fiber) {
					ecmAtBorder.add((ECM) ecmIter); // Add to list of ecm bordering fibers
					fiberNumberTemp.add((Fiber) neighborCheck);
				}
			}
		}
		// Go through the lists of ssc and move to a location that borders muscle border
		for (Object ssc : context) {
			if (ssc instanceof SSC) {
				int index = RandomHelper.nextIntFromTo(0, ecmAtBorder.size() - 1); // draw random number with
																					// randomHelper
				Object ecmRandom = ecmAtBorder.get(index); // get the ecm based on the random number chosen-- index
															// within the list
				((SSC) ssc).setFiberAssoc(((Fiber) (fiberNumberTemp.get(index))).getFiberNumber());
				GridPoint ptECM = grid.getLocation(ecmRandom); // Get the ecm location
				grid.moveTo(ssc, ptECM.getX(), ptECM.getY());
			}
		}
	}

	public static List<Object> getSSCs(Context<Object> context) { // Get a list of all the sscs
		List<Object> sscs = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC) {
				sscs.add(obj);
			}
		}
		return sscs;
	}

	public static List<Object> getActiveSSCs(Context<Object> context) { // get active sscs
		List<Object> sscsActive = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).active >= timeToActive) {
				sscsActive.add(obj);
			}
		}
		return sscsActive;
	}

	public static List<Object> getDiffSSCs(Context<Object> context) { // Get differentiated sscs
		List<Object> sscsDiff = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).differentiated >= 1 && ((SSC) obj).active != 0) { // 2 = myocyte, 3 = fused myocyte on fiber edge 
				sscsDiff.add(obj);
			}
		}
		return sscsDiff;
	}

	public static int getNumActSecretingSSCs(Context<Object> context) { // Get num of secreting SSCs
		List<Object> sscsDiff = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).differentiated >= 1) {
				sscsDiff.add(obj);
			}
		}
		return getActiveSSCs(context).size() - sscsDiff.size(); // returns number of active sscs less differentiated cells on fiber 
	}

	public static int getSenescentSSCs(Context<Object> context) { // Get differentiated sscs
		List<Object> sscsSenescent = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).senescent != 0) { // all differentiated cells
				sscsSenescent.add(obj);
			}
		}
		return sscsSenescent.size();
	}

	public static List<Object> getMyoblastSSCs(Context<Object> context) { // Get myoblasts
		List<Object> sscsMyoblast = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).committed == 2 && ((SSC) obj).differentiated == 0
					&& ((SSC) obj).active != 0) {
				sscsMyoblast.add(obj);
			}
		}
		return sscsMyoblast;
	}

	public static List<Object> getSSCOnFiber(Context<Object> context) {
		List<Object> sscsOnFiber = new ArrayList<Object>();
		for (Object obj : context) {
			if (obj instanceof SSC && ((SSC) obj).onEdge == 1) {
				sscsOnFiber.add(obj);
			}
		}
		return sscsOnFiber;
	}

	public int getSSCOnFiberCount() {
		Context context = ContextUtils.getContext(this);
		return getSSCOnFiber(context).size();
	}

	public int getActiveCount() {
		Context context = ContextUtils.getContext(this);
		return getActiveSSCs(context).size();
	}

	public int getDiffCount() {
		Context context = ContextUtils.getContext(this);
		return getDiffSSCs(context).size();
	}

	public int getMyoblastCount() {
		Context context = ContextUtils.getContext(this);
		return getMyoblastSSCs(context).size();
	}

	public int getPax7Count() {
		Context context = ContextUtils.getContext(this);
		return (getSSCs(context).size() - getDiffCount());
	}

	public void setDivideTime(int divideTime) {
		this.divideTime = divideTime;
	}

	public int getDivideTime() {
		return divideTime;
	}

	public void setDiffTime(int diffTime) {
		this.diffTime = diffTime;
	}

	public int getDiffTime() {
		return diffTime;
	}

	public void setMigrationTime(int migrationTime) {
		this.migrationTime = migrationTime;
	}

	public int getMigrationTime() {
		return migrationTime;
	}

	public void setActive(int active) {
		this.active = active;
	}

	public int getActive() {
		return active;
	}

	public int getDiff() {
		return differentiated;
	}

	public void setOnEdge(int onEdge) {
		this.onEdge = onEdge;
	}

	public int getOnEdge() {
		return onEdge;
	}

	public void setMigrated(int migrated) {
		this.migrated = migrated;
	}

	public int getMigrated() {
		return migrated;
	}

	public void setDaughter(int daughter) {
		this.daughter = daughter;
	}

	public int getDaughter() {
		return daughter;
	}

	public void setCommitted(int committed) {
		this.committed = committed;
	}

	public int getCommitted() {
		return committed;
	}

	public void setSisterAssoc(int sisterAssoc) {
		this.sisterAssoc = sisterAssoc;
	}

	public int getSisterAssoc() {
		return sisterAssoc;
	}

	public void setProteinAdd(int proteinAdd) {
		this.proteinAdd = proteinAdd;
	}

	public int getProteinAdd() {
		return proteinAdd;
	}

	public void setFiberNeedsRep(int fiberNeedsRep) {
		this.fiberNeedsRep = fiberNeedsRep;
	}

	public int getFiberNeedsRep() {
		return fiberNeedsRep;
	}

	public void setMyf5neg(int myf5neg) {
		this.myf5neg = myf5neg;
	}

	public int getMyf5neg() {
		return myf5neg;
	}

	public void setDiff(int differentiated) {
		this.differentiated = differentiated;
	}

	public double getSSCactivation() {
		return sscActivation;
	}

	public double getSSCmigration() {
		return sscMigration;
	}

	public double getSSCproliferation() {
		return sscProliferation;
	}

	public double getSSCquiescence() {
		return sscQuiescence;
	}

	public double getSSCdifferentiation() {
		return sscDifferentiation;
	}

	public void setSenescent(int senescent) {
		this.senescent = senescent;
	}

	public int getSenescent() {
		return senescent;
	}

	public void setFiberAssoc(int fiberAssoc) {
		this.fiberAssoc = fiberAssoc;
	}

	public int getFiberAssoc() {
		return fiberAssoc;
	}

	public double getSSCRecruitSatTemp() {
		return sscRecruitSatTemp;
	}

}
