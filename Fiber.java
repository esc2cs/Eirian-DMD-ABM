/**
 * 
 */
package muscle;

import java.awt.Paint;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;

import bibliothek.gui.dock.station.stack.tab.layouting.Size;
import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.query.space.grid.VNQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.graph.Network;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;

/**
 * @author Kelley Virgilio
 * Defines the class of muscle fibers- located in the grid
 */
public class Fiber {
	
	// FIBER PARAMETERS:
	private Grid<Object> grid;
	private  BufferedGridValueLayer mcpSpatial;
	public static int tick; // keep track of tick count for the active tgf
	private int mfProtein; // amount of protein
	public int fiberNumber;	// fiber number-- all fiber agents in same fiber have same fiber number
	private int border; // 1 == border of fiber/ecm, 0 == not a border fiber elem
	public static int origFiberNumber; // baseline number of fibers
	public static double[] inflamCellsInitial; // inflammatory cell counts
	public int needsRepair; // 0 = no repair needed, 1 = needs repair (if needs growth < -1  + damaged)
	public double needsGrowth; // calculates the z score of the fiber damage (value - mean)/stdev; -1 or less = needs growth
	public int damaged; // = 1 if fiber borders necrosis or low collagen areas-- # represents how many elements of necrosis or low collagen border the fiber
	public int sscCountOnFiber; // number of ssc on fiber
	public static double origCsaMean; // average fiber size at beginning of simulations
	public static double origCsaStDev; // average fiber size standard deviation at beginning of simulations
	private static int fibersRepairing; // toggle to let fibers know they are regrowing
	public static int needRepairCount; // counts the total number of fibers that need repair
	public static double[] origFiberCSA; // orig fiber CSA for each fiber
	public static double[] origFiberNecrosis; // orig amount of necrosis for each fiber
	public static double[] chronicFiberNecrosis; // metric for chronic damage simulations

	// DISEASE STATE PARAMETER
	public static int diseaseState = 0; // 0 = healthy, 1 = young mdx, 2 = adult mdx, 3 = old mdx, 4 = very young mdx
	public static double regenCapacity = 1; // 1 at healthy-- altered at disease states
	public static double asymmSenescent = 0; // 0 at healthy-- altered at disease states
	public static double mdxChronicInflam = 1; //1 at healthy-- altered at disease states
	public static double mdxsSecondNecr = 1; //1 at healthy-- altered at disease states
	public static double fibrobMDX = 1; // 1 at healthy -- altered at disease state
	public static double mdxTnf = 0; // 0 at healthy -- altered at disease state
	public static double mdxIfn = 0; // 0 at healthy -- altered at disease state
	public static double mdxM1mult = 1; // 1 at healthy -- altered at disease state
	public static double mdxM2mult = 1; // 1 at healthy -- altered at disease state
	public static double mdxM1death = 1; // 1 at healthy -- altered at disease state
	public static double mdxM2death = 1; // 1 at healthy -- altered at disease state
	public static double mdxEosinophil = 1; // 1 at healthy -- altered at disease state
	public static double mdxBaseCollagen = 1; // 1 at healthy -- altered at disease state
	public static double macDepletion = 0; // 0 at healthy -- 1 at mac depletion **** MAKE SURE = 0 OTHERWISE
	public static double senescencePa = 0; // 0 at healthy -- 1 at 0% senescence analysis only ** MAKE SURE = 0 OTHERWISE
	public static double necrosisChronic = 1; // set for chronic damage only
	
	// MURPHY 2011 SCALING PARAMETER CHECK
	public static final double pax7Scale = 1;
	public static final double tcf4Scale = 1;

	public Fiber(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, int mfProtein, int fiberNumber, int border, int needsRepair, double needsGrowth, int damaged)
	{
		this.grid = grid;
		this.setMfProtein(mfProtein); // initialize with base level of muscle protein
		this.fiberNumber = fiberNumber; 
		this.border = border;
		this.needsRepair = needsRepair;
		this.needsGrowth = needsGrowth;
		this.damaged = damaged;
		this.mcpSpatial = mcpSpatial;
	}
	
	@ScheduledMethod(start = 2, interval = 1, pick = 1, priority = 1) // Only do this one time for each Class (pick = 1)
	// InflamCell, GrowthFactor class only has one 'agent' in the context
	public void fiberStep(){
		Context context = ContextUtils.getContext(this); // get the context of the fiber
		Neutrophil.neutrophilRecruitment(context);
		Macrophage.macRecruitment(context);
		for (int i = 1; i < origFiberNumber + 1; i++){ // go through each fiber and change the border to red
			List<Object> borderFibers = getFiberBorder(i, context); // get all the borders and set to 1; check if damage surrounds it
			checkFiberSize(context, i); // check if it is big enough
			SSC.sscMigrationChance(context); // chance of SSC migration with every step-- step at every fiber (Scalable with more fibers)
			Fibroblast.fibrobRecruitment(context);// chance of recruiting fibroblasts
			checkFiberTouching(context, i); // check if any fibers are touching and resolve it
		}
		getNumFiberNeedRepair(context); // Get the total number of fibers that need repair
		//System.out.println(tick);
	}
	
	// Initialization
	@ScheduledMethod(start = 1.1, interval = 0, pick = 1) // Complete one time at initialization
	public void fiberInitialStep(){
		Context context = ContextUtils.getContext(this); // get the context of the fiber
		origFiberNumber = getTotalFiberNumber(context); // get total number of fibers at start
		necrosisChronic = 1; // chronic damage only
		fiberCount = 0;
		tick = 0;
		
		// DISEASE PARAMETER SET:// Set disease state parameters
		// HEALTHY
		double necrosisInitial = 35;
		mdxChronicInflam = 1;
		// VERY YOUNG MDX
		if (diseaseState == 4){
			// INITIAL DAMAGE
			necrosisInitial = 35; //
			// STEM CELLS
			asymmSenescent = 3; // 3.5% chance
			// CHRONIC INFLAMMATION
			mdxChronicInflam = 2;
			mdxsSecondNecr = 2; // increases secondary necrosis
			GrowthFactors.mdxTGF = 0; // add TGF-beta in background
//			mdxEosinophil = 1;
			mdxM1mult = 1.5; // 1.5 at very young
			mdxM2mult = .8; // .8 at  very young
			mdxM1death = 1.2; // 1.2 at very young
			mdxM2death = 6; // 6 at very young
			mdxEosinophil = 1.3;
		}
		// YOUNG MDX
		else if (diseaseState == 1){
			// INITIAL DAMAGE
			necrosisInitial = 26;
			// FIBROSIS- COLLAGEN
			mdxBaseCollagen = 1.5;
			List<Object> ecms = ECM.getECM(context);
			for (Object ecm : ecms){
					((ECM) ecm).setCollagen(((ECM) ecm).getCollagen()*mdxBaseCollagen); 
			}
			// STEM CELLS
			// Asymmetric division
			asymmSenescent = 1; // 7% chance
			// CHRONIC INFLAMMATION
			mdxChronicInflam = 1.3; // increases resident macrophages
			mdxsSecondNecr = 1.5; // increases secondary necrosis
			GrowthFactors.mdxTGF = 0; // add TGF-beta in background -- 200
			mdxM2death = .8; // .8 at early  mdx
			// FIBROBLASTS
			fibrobMDX = 1.5;		
		}
		// ADULT MDX
		else if (diseaseState == 2){
			// INITIAL DAMAGE
			necrosisInitial = 29;
			// FIBROSIS- COLLAGEN
			mdxBaseCollagen = 3.0;
			List<Object> ecms = ECM.getECM(context);
			 for (Object ecm : ecms){
			 		((ECM) ecm).setCollagen(((ECM) ecm).getCollagen()*mdxBaseCollagen); 
			 }
			// STEM CELLS
			// Asymmetric division
			asymmSenescent = 2; // 17% chance
			// CHRONIC INFLAMMATION
			mdxChronicInflam = 1.8;
			mdxsSecondNecr = 1.5; // increases secondary necrosis
			GrowthFactors.mdxTGF = 600; // add TGF-beta in background
			mdxM2death = .4; // .4 at adult mdx
			// FIBROBLASTS
			fibrobMDX = 2;			
		}
		// OLD MDX
		else if (diseaseState == 3){
			// INITIAL DAMAGE
			necrosisInitial = 33;
			// FIBROSIS- COLLAGEN
			List<Object> ecms = ECM.getECM(context);
			for (Object ecm : ecms){
				((ECM) ecm).setCollagen(((ECM) ecm).getCollagen()*8); 
			}
			// REGENERATION CAPACITY
			regenCapacity = 3; // increase regenCapacity value
			// STEM CELLS
			// Asymmetric division --> chance of senescence
			asymmSenescent = 1;
			// CHRONIC INFLAMMATION
			mdxChronicInflam = 1.5;
			mdxsSecondNecr = 3; // increases secondary necrosis
			// FIBROBLASTS
			fibrobMDX = 2;
		}

		// get average Fiber CSA and standard deviation at start
		double[] origFiberCSATemp = new double[origFiberNumber];
		double[] origFiberNecrosisTemp = new double[origFiberNumber];
		for (int i = 1; i < origFiberNumber + 1; i++){ // go through each fiber and change the border to red
			getFiberBorder(i, context); // get all the borders and set to 1
			origFiberCSATemp[i-1] = getElemInFiber(i, context).size(); // record the size of each original fiber
		}
		origFiberCSA = origFiberCSATemp;
 		// Calcuate standard deviation and average fiber size for all fibers
	    origCsaMean = StatUtils.mean(origFiberCSATemp, 0, origFiberNumber);
	    origCsaStDev = Math.sqrt(StatUtils.populationVariance(origFiberCSATemp));
	    // chronic damage parameters only:
	    if (InflamCell.chronicDamage == 0){
	    	Necrosis.necrosisInitialize(context,grid, necrosisInitial);
	    }
		// Initialize ssc
		SSC.initialize(origFiberNumber, context, grid); 
		// Initialize fibroblasts	
		Fibroblast.initialize(context, grid, origFiberNumber, fibrobMDX);
		// Initialize inflammatory cells and growth factors
		inflamCellsInitial = InflamCell.initialize(origFiberNumber,  mdxChronicInflam, context);
		double[] inflamCells = getInflamCellInitial(); // Initial values for all the inflammatory cells
	    InflamCell.setInflamCellsIter(inflamCells);
	    GrowthFactors.initializeActiveTgf();
	    // Go back through and set the border- after damage to define which fibers were originally damaged:
	    for (int i = 1; i < origFiberNumber + 1; i++){ // go through each fiber and change the border to red
			getFiberBorder(i, context); // get all the borders and set to 1
			List<Object> elemsInFiber = getElemInFiber(i, context);
			for (Object elems : elemsInFiber){
				if (((Fiber) elems).getDamaged() != 0){ // check if damaged
					origFiberNecrosisTemp[i-1] = 1; // if any of the fibers are marked as damaged- end and go to the next fiber and check
				}
			}
		}
	    origFiberNecrosis = origFiberNecrosisTemp;
	}
	
	
	
	public void getNumFiberNeedRepair(Context<Object> context){ 
		// Get the total number of fibers that need repair
		needRepairCount = 0; // Reset needRepair to 0
		for (int i = 1; i < (origFiberNumber + 1); i++){ // loop through all the fibers
			List<Object> elemInFiber = getElemInFiber(i, context);
			// Get an element in the fiber and check if it needs repair
			if (elemInFiber.size() > 0){
				Object randomElem = elemInFiber.get(0);
				if (((Fiber) randomElem).needsRepair == 1){
					needRepairCount = needRepairCount + 1;
				}
			}
		}
	}
	
	public static double[] getInflamCellInitial(){
		return inflamCellsInitial;
	}
	
	public double addFiberElem(int fiberNumber, Context<Object> context, int removeExt){
		// ADD FIBER ELEMENTS AT MUSCLE REGROWTH
		// IF LOW COLLAGEN  --> REPLACE ECM WITH MUSCLE AND MOVE COLLAGEN TO NEARBY ECM
		// IF HIGH COLLAGEN --> PUSH FIBER AND INCREASE TOTAL CSA OF MUSCLE CROSS SECTION
		// LOW COLLAGEN:
		double lowCollagen = 1.2; // start at 1 and see if any of the values are less
		double collAtRemoveExt = 0.0;
		GridPoint lowCollPt = null;
		Object ecmToRemove = null;
		List<Object> fiberBorderElem = getFiberBorder(fiberNumber, context);
		// Find the neighbor with the lowest collagen, the neighbor ECM cannot ALSO neighbor a fiber with a different fiber number
		List<Object> lowCollAddTemp = new ArrayList<Object>();
		List<Object> allCollAddTemp = new ArrayList<Object>();
		List<Object> adjFiber = new ArrayList<Object>();
		for (Object fiberTemp : fiberBorderElem){ // go through all the fibers at the border and determine which one has an ecm neighbor
				VNQuery<Object> neighborQuery = new VNQuery(grid, fiberTemp, 1, 1); // get 8 potential neighbors
				Iterable<Object> neighborIter = neighborQuery.query();
				// Note: this will always search in a certain order if there are multiple low collagen elements in an area
				for (Object checkNeighbor : neighborIter){ // go through neighbors and find the one with the lowest collagen
					// FIND LOW COLLAGEN NEIGHBOR
					if (checkNeighbor instanceof ECM && ((ECM) checkNeighbor).getCollagen() < lowCollagen){
						// confirm that they do not also border a fiber with a different fiberNumber
						MooreQuery<Object> otherFiberQuery = new MooreQuery(grid, checkNeighbor, 1, 1); // get neighbors
						Iterable<Object> otherFiberIter = otherFiberQuery.query();
						int bordersFiber = 0;
						int bordersEdge = 0;
						for (Object otherFiberCheck : otherFiberIter){ // if the ecm borders a fiber from a different fiber- dont add
							if (otherFiberCheck instanceof Fiber && ((Fiber) otherFiberCheck).getFiberNumber() != ((Fiber) fiberTemp).getFiberNumber()){
								bordersFiber = 1; // go through all neighbors first and check if it borders a fiber
							}
						}
						// check to make sure it is not going to border the exterior grid
						GridPoint ptNeighborCheck2 = grid.getLocation(checkNeighbor);
						if (ptNeighborCheck2.getX() - 1 > 0 && ptNeighborCheck2.getY() - 1 > 0 && ptNeighborCheck2.getX() + 1 < muscleBuilder.gridx && ptNeighborCheck2.getY() + 1 < muscleBuilder.gridy){
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX(),ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX(),ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdge = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()).iterator().hasNext()){
								bordersEdge = 1;
							}
						}
						// Also need to confirm that lowCollagen has a lowCollagen neighbor to put its ECM into ** need TWO adjacent low collagen areas
						if (bordersFiber == 0 && bordersEdge == 0) { // only if the ecm does NOT border another fiber with a different fiber number
							// Check if there is a place to push the low collagen on to
							VNQuery<Object> neighborECMQuery = new VNQuery(grid, ((ECM)checkNeighbor), 1, 1); // get neighbors of this ecm
							Iterable<Object> neighborECMIter = neighborECMQuery.query(); // query the list of agents that are the neighbors
							double lowCollagenTemp2 = 1.2;
							for (Object ecmNeighbor : neighborECMIter){ // find neighbor with lowest ecm and add the ECM here
								if (ecmNeighbor instanceof ECM && ((ECM) ecmNeighbor).getCollagen() < lowCollagenTemp2){
									lowCollAddTemp.add((ECM) checkNeighbor); // List of all possible locations to move with low collagen and not bordering an edge
									lowCollagen = ((ECM) checkNeighbor).getCollagen(); // update low collagen amount and iterate though
								}
							}
						}
					}
					// FIND ECM NEIGHBORS- NOT LOW COLLAGEN-- BUT NEED TO FILL THE HOLES PREFERENTIALLY
					if (checkNeighbor instanceof ECM){
						// confirm that they do not also border a fiber with a different fiberNumber
						MooreQuery<Object> otherFiberQueryAll = new MooreQuery(grid, checkNeighbor, 1, 1); // get neighbors
						Iterable<Object> otherFiberIterAll = otherFiberQueryAll.query();
						int bordersFiberAll = 0;
						int bordersEdgeAll = 0;
						int sameFiberBorder = 0;
						for (Object otherFiberCheck : otherFiberIterAll){ // if the ecm borders a fiber from a different fiber- dont add
							if (otherFiberCheck instanceof Fiber && ((Fiber) otherFiberCheck).getFiberNumber() != ((Fiber) fiberTemp).getFiberNumber()){
								bordersFiberAll = 1; // go through all neighbors first and check if it borders a fiber
							}
							// CHECK TO SEE HOW MANY FIBERS FROM THE SAME FIBER THERE ARE-- FIND THE GREATEST == HOLE
							if (otherFiberCheck instanceof Fiber && ((Fiber) otherFiberCheck).getFiberNumber() == ((Fiber) fiberTemp).getFiberNumber()){
								sameFiberBorder++; // go through all neighbors first and check if it borders a fiber
							}
						}
						((ECM) checkNeighbor).setSameFiberBorder(sameFiberBorder); // set each time- value is not updated unless required
						// check to make sure it is not going to border the exterior grid
						GridPoint ptNeighborCheck2 = grid.getLocation(checkNeighbor);
						if (ptNeighborCheck2.getX() - 1 > 0 && ptNeighborCheck2.getY() - 1 > 0 && ptNeighborCheck2.getX() + 1 < muscleBuilder.gridx && ptNeighborCheck2.getY() + 1 < muscleBuilder.gridy){
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX(),ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()+1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()+1,ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX(),ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()-1).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
							if (!grid.getObjectsAt(ptNeighborCheck2.getX()-1,ptNeighborCheck2.getY()).iterator().hasNext()){
								bordersEdgeAll = 1;
							}
						}
						if (bordersFiberAll == 0 && bordersEdgeAll == 0) { // only if the ecm does NOT border another fiber with a different fiber number
							allCollAddTemp.add((ECM) checkNeighbor); // List of all possible locations to move with low collagen and not bordering an edge
							adjFiber.add((Fiber) fiberTemp); // Also keep track of the adj fiber for the "push out"
						}
					}	
				}
			}
		// After looping through all the fiber edges to find a low collagen neighbor (in a specified fiber number) or the open ECM
		Object ecmToAdd = null;
		if (lowCollAddTemp.size() != 0 && removeExt == 0){
			int randomInt = RandomHelper.nextIntFromTo(0, lowCollAddTemp.size()-1);
			ecmToRemove = lowCollAddTemp.get(randomInt);
			GridPoint newLoc = grid.getLocation(ecmToRemove); // get location
			double collagenAtPt = ((ECM) ecmToRemove).getCollagen();
			// MOVE ECM TO THE NEIGHBOR
			VNQuery<Object> neighborECMQuery = new VNQuery(grid, ecmToRemove, 1, 1); // get neighbors
			Iterable<Object> neighborECMIter = neighborECMQuery.query(); // query the list of agents that are the neighbors
			context.remove(ecmToRemove); // remove ecm
			double collagenNeighbor = 1.2;
			for (Object ecmNeighbor : neighborECMIter){ // find neighbor with lowest ecm and add the ECM here
				if (ecmNeighbor instanceof ECM && ((ECM) ecmNeighbor).getCollagen() < collagenNeighbor){
					ecmToAdd = ecmNeighbor;
					collagenNeighbor = ((ECM) ecmNeighbor).getCollagen(); 
				}
			}
			Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 0, 0, 0, 0); // change to a fiber
			context.add((Fiber) newFiberElem); // change to a fiber
			grid.moveTo(newFiberElem, newLoc.getX(), newLoc.getY());
			if (ecmToAdd != null){ // move collagen to nearby ecm
				((ECM) ecmToAdd).setCollagen(((ECM) ecmToRemove).getCollagen() + collagenNeighbor);
			}
			// NOTE: If ecm to add is null --> It won't move the collagen to a neighbor
		}
		else if (removeExt == 0){
			if (RandomHelper.nextIntFromTo(0,(int)(Math.round(ECM.collagenDensity))) <= 1){ // 1 to density makes chance 100% or 33% for adult; 0 to density changes from 100% to 50% 
			// NO LOW COLLAGEN NEIGHBOR-- NEED TO ADD FIBER AND PUSH OUT ECM
				// INCREASED ENERGY COST TO PUSH OUT-- BASED ON TOTAL COLLAGEN DENSITY
				// NOTE: If removeExtension == 1; then just replace an ecm edge with a fiber instead of pushing out
				Object newDirect = null;
				Object fromDirect = null;
				List<Object> fiberBorderElem2 = getFiberBorder(fiberNumber, context);			
				// Location where fiber will push out
				if (allCollAddTemp.size() > 0){
					int randomInt2 = RandomHelper.nextIntFromTo(0, allCollAddTemp.size()-1);
					newDirect = allCollAddTemp.get(randomInt2);
					fromDirect = adjFiber.get(randomInt2);// get the fiber
					GridPoint newLoc = grid.getLocation(newDirect); // get location of where the fiber is going to push out
					GridPoint fromLoc = grid.getLocation(fromDirect);
					fiberPush(newLoc, fromLoc, fiberNumber); // push out fibers
				}
			}
		}
		else if (removeExt == 1 && allCollAddTemp.size() > 0){
			// Called from ECM.removeExtension--> removeExt == 1
			// Replace the ecm with a fiber
			// Choose the ecm with the highest number of sameFiberBorder
			int sameFiberBorderTemp = 0;
			ecmToRemove = null;
				for (Object ecmAll : allCollAddTemp){
					if (((ECM) ecmAll).getSameFiberBorder() > sameFiberBorderTemp){
						ecmToRemove = ((ECM) ecmAll);
					}
				}
			// Now replace the ECM with the most fiber borders
			GridPoint locationToChange = grid.getLocation(ecmToRemove);
			double collagenAtPt = ((ECM) ecmToRemove).getCollagen();
			context.remove(ecmToRemove); // remove ecm
			collAtRemoveExt = collagenAtPt;
			// Return collagen and it will be replaced where the other ECM was added
			// MOVE ECM TO THE NEIGHBOR
			Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 0, 0, 0, 0); // change to a fiber
			context.add((Fiber) newFiberElem); // change to a fiber
			grid.moveTo(newFiberElem, locationToChange.getX(), locationToChange.getY());
		}
		return collAtRemoveExt;
	}
	
	public void moveByOne(int oldx, int oldy, int newx, int newy){
		while (grid.getObjectAt(oldx,oldy) != null){ // move all the agents
			grid.moveTo(grid.getObjectAt(oldx,oldy), newx, newy);
		}	
	}
	
	public void fiberPush(GridPoint newLoc, GridPoint fromLoc, int fiberNumber){
		// ADD A FIBER ELEMENT AND PUSH THE ENTIRE ROW/COLUMN OUT
		// Calculation direction to grow based on newLoc - fromLoc
		Context context = ContextUtils.getContext(this);
		int xDir = newLoc.getX() - fromLoc.getX();
		int yDir = newLoc.getY() - fromLoc.getY();
		if (xDir != 0){
			// Move in the x direction
			if (xDir == 1){
/////////////////////// POSITIVE X
				// Find the outer ecm edge in this direction
				List<Object> ecmEdge = ECM.defineEdge(context, grid);
				List<Object> ecmEdgeInLine = new ArrayList();
				// reset the ecm edge
				for (Object ecmEdgeTemp : ecmEdge){
					if (grid.getLocation((ECM)ecmEdgeTemp).getY() == newLoc.getY()){
						// find the ecm edge elems in this row
						ecmEdgeInLine.add((ECM)ecmEdgeTemp);
					}
				}
				int tempEdge = 0;
				Object ecmEdgeTemp = null;
				for (Object ecmedges : ecmEdgeInLine){
					// TOP RIGHT --> RIGHT EDGE = HIGH X
					if (grid.getLocation((ECM)ecmedges).getX() > tempEdge){
						tempEdge = grid.getLocation((ECM)ecmedges).getX();
						ecmEdgeTemp = (ECM)ecmedges;
					}
				}
				for (int i = Math.abs(tempEdge - newLoc.getX()) + 1; i > 0; i--){
					// Move the furthest element first
					int oldx = newLoc.getX() + i - 1;
					int oldy = newLoc.getY();
					int newx = newLoc.getX() + i ;
					int newy = oldy;
					if (newx < muscleBuilder.gridx){
						moveByOne(oldx, oldy, newx, newy);
					}
					else {
						return; // Must stop if this is too far at the edge
					}
				}
				// NEED TO GET FIBER NUMBER
				Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 1, 0, 0, 0); // change to a fiber
				context.add((Fiber) newFiberElem); // change to a fiber
				grid.moveTo(newFiberElem, newLoc.getX(), newLoc.getY());
				//System.out.println(newFiberElem);
			}
			else if (xDir == -1){
/////////////////////// NEGATIVE X
				// Find the outer ecm edge in this direction
				List<Object> ecmEdge = ECM.defineEdge(context, grid);
				List<Object> ecmEdgeInLine = new ArrayList();
				// reset the ecm edge
				for (Object ecmEdgeTemp : ecmEdge){
					if (grid.getLocation((ECM)ecmEdgeTemp).getY() == newLoc.getY()){
						// find the ecm edge elems in this row
						ecmEdgeInLine.add((ECM)ecmEdgeTemp);
					}
				}
				int tempEdge = muscleBuilder.gridx;
				Object ecmEdgeTemp = null;
				for (Object ecmedges : ecmEdgeInLine){
					// TOP RIGHT --> RIGHT EDGE = HIGH X
					if (grid.getLocation((ECM)ecmedges).getX() < tempEdge){
						tempEdge = grid.getLocation((ECM)ecmedges).getX();
						ecmEdgeTemp = (ECM)ecmedges;
					}
				}
				for (int i = Math.abs(newLoc.getX() - tempEdge) + 1; i > 0; i--){
					// Move the furthest element first
					int oldx = newLoc.getX() - i + 1;
					int oldy = newLoc.getY();
					int newx = newLoc.getX() - i ;
					int newy = oldy;
					if (newx >= 0){
						moveByOne(oldx, oldy, newx, newy);
					}
					else {
						return; // Must stop if this is too far at the edge
					}
				}
				// NEED TO GET FIBER NUMBER
				Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 1, 0, 0, 0); // change to a fiber
				context.add((Fiber) newFiberElem); // change to a fiber
				grid.moveTo(newFiberElem, newLoc.getX(), newLoc.getY());
				//System.out.println(newFiberElem);
			}
		}
		else if (yDir != 0){
/////////////////////// POSITIVE Y
			if (yDir == 1){
				// Move in positive y
				// Find the outer ecm edge in this direction
				List<Object> ecmEdge = ECM.defineEdge(context, grid);
				List<Object> ecmEdgeInLine = new ArrayList();
				// reset the ecm edge
				for (Object ecmEdgeTemp : ecmEdge){
					if (grid.getLocation((ECM)ecmEdgeTemp).getX() == newLoc.getX()){
						// find the ecm edge elems in this row
						ecmEdgeInLine.add((ECM)ecmEdgeTemp);
					}
				}
				int tempEdge = 0;
				Object ecmEdgeTemp = null;
				for (Object ecmedges : ecmEdgeInLine){
					// TOP RIGHT --> RIGHT EDGE = HIGH X
					if (grid.getLocation((ECM)ecmedges).getY() > tempEdge){
						tempEdge = grid.getLocation((ECM)ecmedges).getY();
						ecmEdgeTemp = (ECM)ecmedges;
					}
				}
				for (int i = Math.abs(tempEdge - newLoc.getY()) + 1; i > 0; i--){
					// Move the furthest element first
					int oldx = newLoc.getX() ;
					int oldy = newLoc.getY()+ i - 1;
					int newx = oldx;
					int newy = newLoc.getY() + i ;
					if (newy <= muscleBuilder.gridy){
						moveByOne(oldx, oldy, newx, newy);
					}
					else {
						return; // Must stop if this is too far at the edge
					}
				}
				// NEED TO GET FIBER NUMBER
				Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 1, 0, 0, 0); // change to a fiber
				context.add((Fiber) newFiberElem); // change to a fiber
				grid.moveTo(newFiberElem, newLoc.getX(), newLoc.getY());
				//System.out.println(newFiberElem);
			}
			else if (yDir == -1){
/////////////////////// NEGATIVE Y
				// Find the outer ecm edge in this direction
				List<Object> ecmEdge = ECM.defineEdge(context, grid);
				List<Object> ecmEdgeInLine = new ArrayList();
				// reset the ecm edge
				for (Object ecmEdgeTemp : ecmEdge){
					if (grid.getLocation((ECM)ecmEdgeTemp).getX() == newLoc.getX()){
						// find the ecm edge elems in this row
						ecmEdgeInLine.add((ECM)ecmEdgeTemp);
					}
				}
				int tempEdge = muscleBuilder.gridy;
				Object ecmEdgeTemp = null;
				for (Object ecmedges : ecmEdgeInLine){
					// TOP RIGHT --> RIGHT EDGE = HIGH X
					if (grid.getLocation((ECM)ecmedges).getY() < tempEdge){
						tempEdge = grid.getLocation((ECM)ecmedges).getY();
						ecmEdgeTemp = (ECM)ecmedges;
					}
				}
				for (int i = Math.abs(newLoc.getY() - tempEdge) + 1; i > 0; i--){
					// Move the furthest element first
					int oldx = newLoc.getX();
					int oldy = newLoc.getY() - i + 1;
					int newx = oldx;
					int newy = newLoc.getY() - i;
					if (newy >= 0 ){
						moveByOne(oldx, oldy, newx, newy);
					}
					else {
						return; // Must stop if this is too far at the edge
					}
				}
				// NEED TO GET FIBER NUMBER
				Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 1, 0, 0, 0); // change to a fiber
				context.add((Fiber) newFiberElem); // change to a fiber
				grid.moveTo(newFiberElem, newLoc.getX(), newLoc.getY());
				//System.out.println(newFiberElem);
			}
		}
	}
	
	public void checkFiberTouching(Context<Object> context, int i){
		// Check if agents in one fiber are touching another fiber
		List<Object> fiberBorderElems = getFiberBorder(i, context); // get all the border fibers
		Object thatFiber = null;
		for (Object thisFiber : fiberBorderElems){ 
			// get the neighbors
			MooreQuery<Object> query = new MooreQuery(grid, thisFiber, 1, 1); // get neighbors
			Iterable<Object> iter = query.query();
			int bordersFiber = 0;
			int bordersEdge = 0;
			for (Object thatFiberTemp : iter){ // if the ecm borders a fiber from a different fiber- dont add
				if (thatFiberTemp instanceof Fiber && ((Fiber) thatFiberTemp).getFiberNumber() != ((Fiber) thisFiber).getFiberNumber()){
					bordersFiber = 1; // go through all neighbors first and check if it borders a fiber
					thatFiber = ((Fiber) thatFiberTemp);
				}
			}
			if (bordersFiber == 1){
				// Change one of them to ECM
				if (RandomHelper.nextIntFromTo(0, 1) < 1){
					GridPoint thisPoint = grid.getLocation((Fiber) thisFiber);
					context.remove((Fiber) thisFiber);
					// Add in ECM
					ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to ECM
					context.add((ECM) newECM);
					grid.moveTo(newECM, thisPoint.getX(), thisPoint.getY()); 
				}
				else{
					GridPoint thatPoint = grid.getLocation((Fiber) thatFiber);
					context.remove((Fiber) thatFiber);
					// Add in ECM
					ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to ECM
					context.add((ECM) newECM);
					grid.moveTo(newECM, thatPoint.getX(), thatPoint.getY()); 
				}
			}
		}

	}
	
	public void checkFiberSize(Context<Object> context, int fiberNumber){
		// At each iteration check if the fiber size is within a stdev of the mean
		if (RandomHelper.nextIntFromTo(0, 5) < 1){ // Only check every couple steps
			List<Object> elemInFibers = getElemInFiber(fiberNumber, context); // get all the elements (agents) in a fiber
			int fiberSize = elemInFibers.size(); // get fiber size
			// Compare to original fiber size of the fiber
			double fiberZ = (fiberSize - origFiberCSA[fiberNumber-1])/origCsaStDev; // compare current fiber to it's initial size
			for (Object fiberTemp : elemInFibers){
				((Fiber) fiberTemp).setNeedsGrowth(fiberZ);
				// getDamaged = 1 signals that the muscle was originally injured
				// NORMAL SINGLE BOUT OF NECROSIS
				if (InflamCell.chronicDamage == 0){
					if (((Fiber) fiberTemp).getNeedsGrowth() < .1  && origFiberNecrosis[fiberNumber-1] != 0){ // if the fiber is smaller than its original size -->  needs repair
						// only for macrophage depletion! --> turn off signal if filled with fibrosis
						if (Fiber.macDepletion == 1 && ECM.collagenDensity > 1.5 && RandomHelper.nextIntFromTo(0, 1) < 1){
							((Fiber) fiberTemp).setNeedsRepair(0);
						}
						else {
							((Fiber) fiberTemp).setNeedsRepair(1);
						} // macrophage depletion test
					}
					else {
						((Fiber) fiberTemp).setNeedsRepair(0);
					}
				}
				// CHRONIC, REPETITIVE DAMAGE ONLY
				if (InflamCell.chronicDamage == 1){
					//System.out.println(chronicFiberNecrosis.length);
					//System.out.println(fiberNumber);

					if (((Fiber) fiberTemp).getNeedsGrowth() < .1  && chronicFiberNecrosis[fiberNumber-1] != 0){ // if the fiber is smaller than its original size -->  needs repair
						// only for macrophage depletion! --> turn off signal if filled with fibrosis
						if (Fiber.macDepletion == 1 && ECM.collagenDensity > 1.5 && RandomHelper.nextIntFromTo(0, 1) < 1){
							((Fiber) fiberTemp).setNeedsRepair(0);
						}
						else {
							((Fiber) fiberTemp).setNeedsRepair(1);
						}
					}
					else {
						((Fiber) fiberTemp).setNeedsRepair(0);
					}
				}
			}
		}
	}

	
	public static List<Object> getFiberElems(Context<Object> context){
		// Get all the fiber elements (agents)
		List<Object> fibers = new ArrayList<Object>(); // create a list of all the fiber agents
		if (context != null){
			for (Object obj : context){ 
				if (obj instanceof Fiber){
					fibers.add(obj);
				}
			} 
			return fibers;
		}
		else return null;
	}
	
	public static List<Object> getElemInFiber(int fiberNumber, Context<Object> context){
		// Create an individual list of all the elements in a specific fiber	
		List<Object> elemInFiber = new ArrayList<Object>();
		List<Object> fibers = getFiberElems(context); // get all the fibers in the context
		for (Object fiberTemp : fibers ){ // go through them all and find the max value
			if (((Fiber) fiberTemp).getFiberNumber() == fiberNumber){
				elemInFiber.add((Fiber) fiberTemp);
			}
		}
		return elemInFiber;
	}
	
	private static int fiberCount = 0;
	@ScheduledMethod(start = 1, interval = 0) // Completed at initialization only
	public void fiberGroup(){ 	// Assign variable for fiber number- User recursive method to solve
		// Iterates through all the fibers in the context and assigns a fiber number
		if (this.getFiberNumber() != 0)
			return;
		// This is where the fiber group counting starts from, recursive method will solve it
		fiberCount = fiberCount + 1;
		Context context = ContextUtils.getContext(this); // get the context of the fiber
		MooreQuery<Object> query = new MooreQuery(grid, this, 1, 1); // Moore finds the 8 neighbors within the extent of the area defined
		Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
		List<Object> fiberNeighbors = new ArrayList<Object>(8); // create a list of neighbors to pass to recursive
		for(Object neighbors : iter){
			if (neighbors instanceof Fiber){
				((Fiber) neighbors).setFiberNumber(fiberCount);
				fiberNeighbors.add((Fiber) neighbors); // create a list of the fiberNeighbors
			}
		}
		if (fiberNeighbors.size() != 0){  
			fiberRecursiveGroup(context, fiberNeighbors, fiberCount); // if there are neighbors that are fibers- send to recursive function
		}
	}
	
	public void fiberRecursiveGroup(Context<Object> context, List<Object> fiberNeighbors, int fiberCount){
		List<Object> fiberNeighborsNew = new ArrayList<Object>(); // create a list of neighbors of neighbors
		for (Object neighbors : fiberNeighbors){
			// Go through the list of the neigbors and get their neighbors
			MooreQuery<Object> query = new MooreQuery(grid, (Fiber) neighbors, 1, 1); // Moore finds the 8 neighbors within the extent of the area defined
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			for(Object neighborsNew : iter){
				if (neighborsNew instanceof Fiber && ((Fiber) neighborsNew).getFiberNumber() == 0){
					((Fiber) neighborsNew).setFiberNumber(fiberCount);
					fiberNeighborsNew.add((Fiber) neighborsNew); // create a list of the fiberNeighbors- only want to add if it is a fiber
				}
			}
		}
		if (fiberNeighborsNew.size() == 0)
			return; // return if there are no neighbors with fibers, else continue
		fiberNeighbors = fiberNeighborsNew;
		fiberRecursiveGroup(context, fiberNeighbors, fiberCount);
	}
	
	public static int getTotalFiberNumber(Context<Object> context){ 
		// Get the max fiberNumber
		List<Object> fibers = getFiberElems(context); // get all the fibers in the context
		int totalFiberNumber = 0;
		for (Object fiberTemp : fibers ){ // go through them all and find the max value
			if (totalFiberNumber < ((Fiber) fiberTemp).getFiberNumber()){
				totalFiberNumber = ((Fiber) fiberTemp).getFiberNumber();
			}
		}
		return totalFiberNumber;
	}
	
	public List<Object> getFiberBorder(int fiberNumber, Context<Object> context){
		// Get the boundary of the fibers for a specific fiber number
		List<Object> elemInFiber = getElemInFiber(fiberNumber, context); // iterate through and get the elem in each fiber number
		List<Object> fiberBorderElem = new ArrayList<Object>();
		int sscCountTemp = 0;
		int getDamagedTemp = 0;
		for (Object elems : elemInFiber){ // go through each element within a specific fiber number
			((Fiber) elems).setBorder(0); // reset to 0
			if (((Fiber) elems).getDamaged() != 0){
				getDamagedTemp = getDamagedTemp + 1; // shows that at least one of the elems in the fiber is damaged
			}
			// Check if there are SSCs on the fiber and set SSC count number for all the elements in that fiber
			if (((Fiber) elems).getsscCountOnFiber() > sscCountTemp){ // find the max sscCounton Fiber and reset for all fibers in the element
				sscCountTemp = ((Fiber) elems).getsscCountOnFiber();
			}
			VNQuery<Object> query = new VNQuery(grid, (Fiber) elems, 1, 1); // Find the 4 neighbors and determine if the fiber borders ecm
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			for(Object neighborCheck : iter){
				if (neighborCheck instanceof ECM || neighborCheck instanceof Necrosis ){ // check if it borders ecm or necrosis
					fiberBorderElem.add((Fiber) elems); // if the neighbor is ECM then add the original fiber to the border group
					((Fiber) elems).setBorder(1); // set border to 1; 1 = border, 0 = no border
				}
			}
		}
		for(Object fiberElems : elemInFiber){ // set needs repair to 1 for all the elements
			((Fiber) fiberElems).setsscCountOnFiber(sscCountTemp);
			if (getDamagedTemp > 0){
				((Fiber) fiberElems).setDamaged(1);
			}
		}
		return fiberBorderElem;
	}
	
	public static int getOrigFiberNumber(){
		return origFiberNumber;
	}
	
	public int getMfProtein() {
		return mfProtein;
	}
	public void setMfProtein(int mfProtein) {
		this.mfProtein = mfProtein;
	}
	
	public void setFiberNumber(int fiberNumber){
		this.fiberNumber = fiberNumber;
	}
	
	public void setNeedsRepair(int needsRepair){
		this.needsRepair = needsRepair;
	}
	
	public int getNeedsRepair(){
		return needsRepair;
	}
	
	public void setNeedsGrowth(double needsGrowth){
		this.needsGrowth = needsGrowth;
	}
	
	public double getNeedsGrowth(){
		return needsGrowth;
	}
	public void setDamaged(int damaged){
		this.damaged = damaged;
	}
	
	public int getDamaged(){
		return damaged;
	}
	
	public static void setFibersRepairing(int fibersRepairingOn){
		fibersRepairing = fibersRepairingOn;
	}
	
	public static int getFibersRepairing(){
		return fibersRepairing;
	}
	
	public void setsscCountOnFiber(int sscCountOnFiber){
		this.sscCountOnFiber = sscCountOnFiber;
	}
	
	public int getsscCountOnFiber(){
		return sscCountOnFiber;
	}
	
	public int getFiberNumber(){
		return fiberNumber;
	}
	
	public void setBorder(int border){
		this.border = border;
	}
	
	public int getBorder(){
		return border;
	}
	
}
