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
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.util.collections.IndexedIterable;
import repast.simphony.util.collections.ListIndexedIterable;
import repast.simphony.util.collections.RandomIterable;
import repast.simphony.valueLayer.BufferedGridValueLayer;

/**
 * @ Kelley Virgilio
 * Necrosis class defines the areas that are necrotic. Replaces the fiber or ecm patch
 *
 */
public class Necrosis {
	
	private Grid<Object> grid;
	private static BufferedGridValueLayer mcpSpatial;
	private int secondary; // defines secondary necrosis
	private int age; // tracks age of necrosis- counter
	public static double initialBurst;
	
	public Necrosis(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, int secondary, int age)
	{
		this.secondary = secondary;
		this.grid = grid;
		this.age = age;
		this.mcpSpatial = mcpSpatial;
	}
	

	@ScheduledMethod(start = 2, interval = 1, priority = 1) 
	public void necrosisRemove(){
		Context context = ContextUtils.getContext(this);
		this.setAge(this.getAge() + 1); // age the necrosis
		initialBurst = getInitialBurstNecrotic(context);
	}

	public static void necrosisBehaviors(Context<Object> context, Grid<Object> grid, double[] inflamCells, int totalFiberNumber, double[] growthFactors){		
		// inflammCells: 0 RM; 1 N; 2 Na; 3 M1; 4 M1ae; 5 M1de; 6 M2
		// Neutrophils and M1-debris eating get rid of necrosis but release ROS --> secondary necrosis	
		/*
		 * 
		 * commented out from here to 107
		 */
		double neutrophils = Math.ceil(inflamCells[1])*50; // round neutrophils to nearest whole number
		double m1de = Math.ceil(Math.ceil(inflamCells[5]) + GrowthFactors.m1MacAdded*.3)*50; // round neutrophils to nearest whole number
		// double m1de = 0;// Math.ceil(inflamCells[5])*50;
		// M1 DEBRIS EATING MACROPHAGES- REMOVE NECROSIS
		  for (int i = 0; i < m1de; i++){
			  // *** CHECK IF THERE IS ANY NECROSIS IN THE MUSCLE
			List<Object> necrosisCheck = getNecrosis(context);
			if (necrosisCheck.size() > 0){  
			  if (RandomHelper.nextIntFromTo(0, 1) < 2){ // Each m1de has a 1 in 10 chance of removing necrosis ***changed to remove always 
				// REMOVE NECROSIS
				List<Necrosis> necrosisAtEdge = new ArrayList<Necrosis>();
				List<Object> necrosis = getNecrosis(context); // Get all the necrosis elements at each step
				// Make a list of all the necrosis that borders ECM
				for (Object necrosisIter : necrosis){ // For all the necrosis elements- get the neighbors
					VNQuery<Object> query = new VNQuery(grid, necrosisIter, 1, 1); // Find the 4 neighbors and determine if the necrosis borders ecm
					Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
					for (Object neighborCheck : iter){ // go through the list of neighbors
						if (neighborCheck instanceof ECM){
							necrosisAtEdge.add((Necrosis) necrosisIter); // Add to list of necrosis elements bordering ECM
						}
					}
				}
				// Now choose a random necrosis element from the list of necrosisAtEdge to change to ECM
				// CHECK IF THERE IS NECROSIS AT EDGE: ELSE CHOOSE ANY NECROSIS		
				if (necrosisAtEdge.size() > 0){
					int index = RandomHelper.nextIntFromTo(0,  necrosisAtEdge.size() - 1);
					Necrosis necrosisRandom = necrosisAtEdge.get(index);
					ECM newECM = new ECM(mcpSpatial, grid, 0.1, 0, 0); // Create a new ECM element currently just add small amount of collagen
					//ECM newECM = new ECM(grid, .1, 0, 0); // add a lot of collagen-- removed by MMPs
					GridPoint pt = grid.getLocation(necrosisRandom); // get location
					context.remove(necrosisRandom); // remove the necrosis and change to ECM
					context.add(newECM); // Change to ecm for now
					grid.moveTo(newECM, pt.getX(), pt.getY()); // Move ECM to wear the necrosis was removed
				}
				else {
					List<Object> necrosisRemain = getNecrosis(context); // get the remaining necrosis
					int index = RandomHelper.nextIntFromTo(0,  necrosisRemain.size() - 1);
					Object necrosisRandomTemp = necrosisRemain.get(index);
					ECM newECM = new ECM(mcpSpatial, grid, 0.1, 0, 0); // Create a new ECM element currently just add small amount of collagen
					GridPoint pt = grid.getLocation(necrosisRandomTemp); // get location
					context.remove(necrosisRandomTemp); // remove the necrosis and change to ECM
					context.add(newECM); // Change to ecm for now
					grid.moveTo(newECM, pt.getX(), pt.getY()); // Move ECM to wear the necrosis was removed
				}
			  }
			}
		} //*/
		
		// SECONDARY NECROSIS- based on ROS
		double ros = growthFactors[28];
		// Note: May need to adjust this parameter as ROS values change with increased cell counts
		//all of this was commented out
		if (ros > 0.001*totalFiberNumber && 1 > RandomHelper.nextIntFromTo(0,10)){ // ROS + chance 
			for (int i = 0; i < ros*Fiber.mdxsSecondNecr; i++){
			// Add necrosis next to a location of necrosis-- change color of secondary necrosis to orange
			List<Object> fibers = Fiber.getFiberElems(context); // get the fibers
			List<Fiber> fiberNearNecrosis = new ArrayList<Fiber>();
			for (Object fiberIter : fibers){
				VNQuery<Object> query = new VNQuery(grid, fiberIter, 1, 1); // Find the 4 neighbors and determine if the necrosis borders ecm
				Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
				for (Object neighborCheck : iter){ // go through the list of neighbors
					if (neighborCheck instanceof Necrosis){
						fiberNearNecrosis.add((Fiber) fiberIter); // Add to list of fibers boarding necrosis
					}
				}
			}
			// Choose a random fiber near necrosis and change it to secondary necrosis
			if (fiberNearNecrosis.size() > 0){
			int index = RandomHelper.nextIntFromTo(0,  fiberNearNecrosis.size() - 1);
			Fiber fiberRandom = fiberNearNecrosis.get(index);
			Necrosis newNecrosis = new Necrosis(mcpSpatial, grid, 1, 0); // Create a new Necrotic element, marked 'secondary'
			GridPoint pt = grid.getLocation(fiberRandom); // get location
			context.remove(fiberRandom); // remove the fiber and change to secondary necrosis
			context.add(newNecrosis); // Add new necrosis
			grid.moveTo(newNecrosis, pt.getX(), pt.getY()); 
			}
			}
		}
	}
	
	//@ScheduledMethod(start = 2, interval = 24, priority = 1) 
	// Only for chronic damage simulations
	public static void chronicDamage(Context<Object> context, Grid<Object> grid, double chronicAmount){
		// Apply a low level chronic damage
		double necrosisInitial = chronicAmount/100.; // parameter held as the percent
		List<Object> fibers = new ArrayList<Object>();
		fibers = Fiber.getFiberElems(context);
		double initialNecrFibtemp = fibers.size()*necrosisInitial;
		double initialNecrFib = Math.floor(initialNecrFibtemp); // TOTAL NUMBER OF FIBERS THAT SHOULD BE NECROTIC AT START
		double currentNecr = 0;
		while (currentNecr < initialNecrFib){ // MAKE FIBERS NECROTIC UNTIL IT HAS REACHED THE INPUT VALUE OF NECROSIS
			fibers = Fiber.getFiberElems(context); // RESET	FIBER LIST AT EACH ITERATION FOR REMAINING FIBERS
			int index = RandomHelper.nextIntFromTo(0,  fibers.size() - 1); // DRAW RANDOM NUMBER FROM THE NUMBER OF FIBERS
			Object fiberRandom = fibers.get(index); // get a random border fiber to start with
			while (((Fiber) fiberRandom).getBorder() == 0){
				int newIndex = RandomHelper.nextIntFromTo(0,  fibers.size() - 1); // get a new random number if not a border fiber
				Object fiberRandomTemp = fibers.get(newIndex); // get a random border fiber to start with
				fiberRandom = fiberRandomTemp;
			}
			// CHANGE FIBERS TO NECROTIC
			Necrosis necrosis = new Necrosis(mcpSpatial, grid, 0, 0);
			GridPoint ptfiberRandom = grid.getLocation(fiberRandom); // GET THE ORIGINAL FIBER LOCATION			
			context.remove(fiberRandom); // REMOVE FIBER FROM THE CONTEXT
			context.add(necrosis); // ADD THE NECROSIS
			grid.moveTo(necrosis, ptfiberRandom.getX(), ptfiberRandom.getY()); // MOVE INTO CONTEXT IN PLACE OF OLD FIBER
			currentNecr = currentNecr + 1; // ADD ONE TO CURRENT NECROTIC
			// EXPAND NECROSIS OUT FROM NECROTIC AREA
			/*
			for (int i = 0; i < RandomHelper.nextIntFromTo(0,  2000); i++){ // TO DO: reset max to be dependent on the size of the fibers (once defined)
				fibers = Fiber.getFiberElems(context); // GET THE NEW LIST OF FIBERS
				MooreQuery<Object> query = new MooreQuery(grid, necrosis, 2, 2); // Von Neumann query finds the 4 neighbors; Moore finds the 8 neighbors within the extent of the area defined
				Iterable<Object> iter = query.query();
				int randomNumber = RandomHelper.nextIntFromTo(0,  2);
				int count = 0;
				// FIX CONCURRENT MOD
				// Create a temp list and remove after
				List<Object> fiberNeighborsTemp = new ArrayList<Object>();
				for (Object fiberNeighbors : iter){
					if (fiberNeighbors instanceof Fiber){
						fiberNeighborsTemp.add((Fiber) fiberNeighbors); // creates a new list of fiberNeighbors
					}
				}
				for (Object fiberNeighborChange : fiberNeighborsTemp){
					Necrosis necrosisNeighbor = new Necrosis(mcpSpatial, grid, 0, 0);
					GridPoint ptNeighbor = grid.getLocation(fiberNeighborChange); // GET THE ORIGINAL FIBER LOCATION			
					context.remove(fiberNeighborChange); // REMOVE FIBER FROM THE CONTEXT
					context.add(necrosisNeighbor); // ADD THE NECROSIS
					grid.moveTo(necrosisNeighbor, ptNeighbor.getX(), ptNeighbor.getY()); // MOVE INTO CONTEXT IN PLACE OF OLD FIBER
					currentNecr = currentNecr + 1;
					// Choose a random neighbor from Fiber so it builds out in a random fashion
					if (count == randomNumber){
						necrosis = necrosisNeighbor; // reset to one of the fibers and get it's neighbors
					}
					count = count + 1;
				}
		} */
			
			// new code to replace lines 159-186
			for (int i = 0; i < RandomHelper.nextIntFromTo(0, 1); i++){ // this line tells how large the cluster of damage will be
			fibers = Fiber.getFiberElems(context); // GET THE NEW LIST OF FIBERS
			MooreQuery<Object> query = new MooreQuery(grid, necrosis, 2, 2); // Von Neumann query finds the 4 neighbors; Moore finds the 8 neighbors within the extent of the area defined
			Iterable<Object> iter = query.query();
			int randomNumber = RandomHelper.nextIntFromTo(0,  2);
			int count = 0;
			// FIX CONCURRENT MOD
			// Create a temp list and remove after
			List<Object> fiberNeighborsTemp = new ArrayList<Object>();
			for (Object fiberNeighbors : iter){
				if (fiberNeighbors instanceof Fiber){
					fiberNeighborsTemp.add((Fiber) fiberNeighbors); // creates a new list of fiberNeighbors
				}
			}
			for (Object fiberNeighborChange : fiberNeighborsTemp){
				Necrosis necrosisNeighbor = new Necrosis(mcpSpatial, grid, 0, 0);
				GridPoint ptNeighbor = grid.getLocation(fiberNeighborChange); // GET THE ORIGINAL FIBER LOCATION			
				context.remove(fiberNeighborChange); // REMOVE FIBER FROM THE CONTEXT
				context.add(necrosisNeighbor); // ADD THE NECROSIS
				grid.moveTo(necrosisNeighbor, ptNeighbor.getX(), ptNeighbor.getY()); // MOVE INTO CONTEXT IN PLACE OF OLD FIBER
				currentNecr = currentNecr + 1;
				// Choose a random neighbor from Fiber so it builds out in a random fashion
				if (count == randomNumber){
					necrosis = necrosisNeighbor; // reset to one of the fibers and get it's neighbors
				}
				count = count + 1;
			}
	} 
			
	
			
	 }
		// After applying necrosis- determine which fibers are considered damaged:
		// Iterate through each fiber and check if it is damaged:
		for (int i = 1; i < Fiber.origFiberNumber + 1; i++){
			List<Object> elemInFiber = Fiber.getElemInFiber(i, context); // iterate through and get the elem in each fiber number
			List<Object> fiberBorderElem = new ArrayList<Object>();
			int fiberDamaged = 0;
			for (Object elems : elemInFiber){ // go through each element within a specific fiber number
				VNQuery<Object> query = new VNQuery(grid, (Fiber) elems, 1, 1); // Find the 4 neighbors and determine if the fiber borders ecm
				Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
				for(Object neighborCheck : iter){
					if ((neighborCheck instanceof ECM && ((ECM) neighborCheck).getCollagen() < 1) || neighborCheck instanceof Necrosis ){
						// If the edge borders necrosis, or past necrosis with low collagen-- mark this is a fiber that needs repair
						fiberDamaged = fiberDamaged + 1;
					}
				}
			}
			if (fiberDamaged > 1){
					for(Object fiberElems : elemInFiber){ // set needs repair to 1 for all the elements
						((Fiber) fiberElems).setDamaged(1); // set damaged to 1 -- never reset to 0
					}
			}
		}
	}
	
	public static void necrosisInitialize(Context<Object> context, Grid<Object> grid, double necroticAmount){
		// INITIALIZE WITH A SET AMOUNT OF DAMAGE IN THE MUSCLE- BASED ON INPUT PARAMETER
		double necrosisInitial = necroticAmount/100.; // parameter held as the percent
		List<Object> fibers = new ArrayList<Object>();
		fibers = Fiber.getFiberElems(context);
		List<Object> ecms = new ArrayList<Object>();
		ecms = ECM.getECM(context);
		double initialNecrFibtemp = (fibers.size())*necrosisInitial; // % of fiber only--- better for simulations for force loss == % fiber (assumption)
		double initialNecrFib = Math.floor(initialNecrFibtemp); // TOTAL NUMBER OF FIBERS THAT SHOULD BE NECROTIC AT START
		double currentNecr = 0;
		while (currentNecr < initialNecrFib){ // MAKE FIBERS NECROTIC UNTIL IT HAS REACHED THE INPUT VALUE OF NECROSIS
			fibers = Fiber.getFiberElems(context); // RESET	FIBER LIST AT EACH ITERATION FOR REMAINING FIBERS
			int index = RandomHelper.nextIntFromTo(0,  fibers.size() - 1); // DRAW RANDOM NUMBER FROM THE NUMBER OF FIBERS
			Object fiberRandom = fibers.get(index); // get a random border fiber to start with
			while (((Fiber) fiberRandom).getBorder() == 0){
				int newIndex = RandomHelper.nextIntFromTo(0,  fibers.size() - 1); // get a new random number if not a border fiber
				Object fiberRandomTemp = fibers.get(newIndex); // get a random border fiber to start with
				fiberRandom = fiberRandomTemp;
			}
			// CHANGE FIBERS TO NECROTIC
			Necrosis necrosis = new Necrosis(mcpSpatial, grid, 0, 0);
			GridPoint ptfiberRandom = grid.getLocation(fiberRandom); // GET THE ORIGINAL FIBER LOCATION			
			context.remove(fiberRandom); // REMOVE FIBER FROM THE CONTEXT
			context.add(necrosis); // ADD THE NECROSIS
			grid.moveTo(necrosis, ptfiberRandom.getX(), ptfiberRandom.getY()); // MOVE INTO CONTEXT IN PLACE OF OLD FIBER
			currentNecr = currentNecr + 1; // ADD ONE TO CURRENT NECROTIC
			// EXPAND NECROSIS OUT FROM NECROTIC AREA
			for (int i = 1; i < RandomHelper.nextIntFromTo(0,  2000); i++){ // TO DO: reset max to be dependent on the size of the fibers (once defined)
				fibers = Fiber.getFiberElems(context); // GET THE NEW LIST OF FIBERS
				MooreQuery<Object> query = new MooreQuery(grid, necrosis, 1, 1); // Von Neumann query finds the 4 neighbors; Moore finds the 8 neighbors within the extent of the area defined
				// *** changed from 2,2 area to 1,1
				Iterable<Object> iter = query.query();
				int randomNumber = RandomHelper.nextIntFromTo(0, 2);
				int count = 0;
					for (Object neighbors : iter){
						if (neighbors instanceof Fiber){
							Necrosis necrosisNeighbor = new Necrosis(mcpSpatial, grid, 0, 0);
							GridPoint ptNeighbor = grid.getLocation(neighbors); // GET THE ORIGINAL FIBER LOCATION			
							context.remove(neighbors); // REMOVE FIBER FROM THE CONTEXT
							context.add(necrosisNeighbor); // ADD THE NECROSIS
							grid.moveTo(necrosisNeighbor, ptNeighbor.getX(), ptNeighbor.getY()); // MOVE INTO CONTEXT IN PLACE OF OLD FIBER
							currentNecr = currentNecr + 1;
							// Choose a random neighbor from Fiber so it builds out in a random fashion
							if (count == randomNumber){
								necrosis = necrosisNeighbor; // reset to one of the fibers and get it's neighbors
							}
							count = count + 1;
						}
					}
		}
	 }
		// After applying necrosis- determine which fibers are considered damaged:
		// Iterate through each fiber and check if it is damaged:
		for (int i = 1; i < Fiber.origFiberNumber + 1; i++){
			List<Object> elemInFiber = Fiber.getElemInFiber(i, context); // iterate through and get the elem in each fiber number
			List<Object> fiberBorderElem = new ArrayList<Object>();
			int fiberDamaged = 0;
			for (Object elems : elemInFiber){ // go through each element within a specific fiber number
				VNQuery<Object> query = new VNQuery(grid, (Fiber) elems, 1, 1); // Find the 4 neighbors and determine if the fiber borders ecm
				Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
				for(Object neighborCheck : iter){
					if ((neighborCheck instanceof ECM && ((ECM) neighborCheck).getCollagen() < 1) || neighborCheck instanceof Necrosis ){
						// If the edge borders necrosis, or past necrosis with low collagen-- mark this is a fiber that needs repair
						fiberDamaged = fiberDamaged + 1;
					}
				}
			}
			if (fiberDamaged > 1){
					for(Object fiberElems : elemInFiber){ // set needs repair to 1 for all the elements
						((Fiber) fiberElems).setDamaged(1); // set damaged to 1
					}
			}
		}
	}
	
	// Get rid of necrotic element and replace with ECM
	public static void eatNecrosis(Context<Object> context, Grid<Object> grid, Object necrosisToEat) {
		ECM newECM = new ECM(mcpSpatial, grid, 0.1, 0, 0); // Create a new ECM element currently just add small amount of collagen
		GridPoint pt = grid.getLocation(necrosisToEat); // get location
		context.remove(necrosisToEat); // remove the necrosis and change to ECM
		context.add(newECM); // Change to ecm for now
		grid.moveTo(newECM, pt.getX(), pt.getY()); // Move ECM to where the necrosis was removed
	}
	
	public static List<Object> getNecrosis(Context<Object> context){ // Get a list of all the fibroblasts
		List<Object> necrosis = new ArrayList<Object>(); // create a list of all the fibroblast agents
		if (context != null){
			for (Object obj : context){ 
				if (obj instanceof Necrosis){
					necrosis.add(obj);
				}
			}
		} 
		return necrosis;
	}
	
	public static double getPercentNecrotic(Context<Object> context){
		List<Object> necrosis =  getNecrosis(context);// Get amount of necrosis
		List<Object> fibers = Fiber.getFiberElems(context);
		List<Object> ecms = ECM.getECM(context);
		double necrosisSize = necrosis.size();
		double ecmSize = ecms.size();
		double fiberSize = fibers.size();
		return necrosisSize/(necrosisSize + fiberSize + ecmSize); // percent necrotic muscle 
	}
	
	public static double getRecentPercentNecrotic(Context<Object> context){
		// Calculates the recent percent necrotic --> amount in last 72 hours
		int necroticSizeTemp = 0; // tracks the amount of necrotic fibers necrosed in last 72 hours
		List<Object> necrosis =  getNecrosis(context);// Get amount of necrosis
		for (Object necroticRec : necrosis){
			if (((Necrosis) necroticRec).getAge() < 72){
				necroticSizeTemp++ ;
			}
		}
		List<Object> fibers = Fiber.getFiberElems(context);
		List<Object> ecms = ECM.getECM(context);
		double ecmSize = ecms.size();
		double fiberSize = fibers.size();
		return necroticSizeTemp/(necroticSizeTemp + fiberSize + ecmSize); // percent necrotic muscle 
	}
	
	public static double getInitialBurstNecrotic(Context<Object> context){
		// Calculates the initial percent necrotic --> amount in first 24 hours- after which it does not secrete hgf
		int necroticSizeTemp = 0; // tracks the amount of necrotic fibers necrosed in first 24 hours
		List<Object> necrosis =  getNecrosis(context);// Get amount of necrosis
		for (Object necroticRec : necrosis){
			if (((Necrosis) necroticRec).getAge() == 2) {// && ((Necrosis) necroticRec).getSecondary() == 0){
				necroticSizeTemp++ ;
			}
		}
		List<Object> fibers = Fiber.getFiberElems(context);
		List<Object> ecms = ECM.getECM(context);
		double ecmSize = ecms.size();
		double fiberSize = fibers.size();
		return necroticSizeTemp/(necroticSizeTemp + fiberSize + ecmSize); // percent necrotic muscle 
	}
	
	public static List<Necrosis> getNewNecrosis(Context<Object> context){
		// Calculates the initial percent necrotic --> amount in first 24 hours- after which it does not secrete hgf
		int necroticSizeTemp = 0; // tracks the amount of necrotic fibers necrosed in first 24 hours
		List<Object> necrosis =  getNecrosis(context);// Get amount of necrosis
		List<Necrosis> newNecrosis = new ArrayList<Necrosis>();
		for (Object necroticRec : necrosis){
			if (((Necrosis) necroticRec).getAge() < 2){
				newNecrosis.add((Necrosis) necroticRec);
			}
		}
		return newNecrosis; // percent necrotic muscle 
	}
	
	public static List<Necrosis> getNecrosisAtEdge(Context<Object> context, Grid<Object> grid) {
		List<Necrosis> necrosisAtEdge = new ArrayList<Necrosis>();
		List<Object> necrosis = getNecrosis(context); // Get all the necrosis elements at each step
		// Make a list of all the necrosis that borders ECM
		for (Object necrosisIter : necrosis){ // For all the necrosis elements- get the neighbors
			VNQuery<Object> query = new VNQuery(grid, necrosisIter, 1, 1); // Find the 4 neighbors and determine if the necrosis borders ecm
			Iterable<Object> iter = query.query(); // query the list of agents that are the neighbors
			for (Object neighborCheck : iter){ // go through the list of neighbors
				if (neighborCheck instanceof ECM){
					necrosisAtEdge.add((Necrosis) necrosisIter); // Add to list of necrosis elements bordering ECM
				}
			}
		}
		return necrosisAtEdge;
	}
	
	public void setSecondary(int secondary){
		this.secondary = secondary;
	}
	
	public int getSecondary(){
		return secondary;
	}
	public void setAge(int age){
		this.age = age;
	}
	
	public int getAge(){
		return age;
	}
}

	

