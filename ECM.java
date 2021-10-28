/**
 * 
 */
package muscle;

import java.util.ArrayList;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.query.space.grid.VNQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.BufferedGridValueLayer;

/**
 * @author Kelley Virgilio
 *
 */
public class ECM {

	private Grid<Object> grid;
	private BufferedGridValueLayer mcpSpatial;
	private double collagen;
	private static ContinuousSpace<Object> space;
	private int ecmEdge;
	private int sameFiberBorder;
	public static double collagenDensity;

	public ECM(BufferedGridValueLayer mcpSpatial, Grid<Object> grid, double collagen, int ecmEdge, int sameFiberBorder) {
		this.grid = grid;
		this.collagen = collagen; // initialize with base level of collagen
		this.ecmEdge = ecmEdge; // marks the outer edge of the ecm
		this.sameFiberBorder = sameFiberBorder; // records how many fibers the ecm borders --> greater number ==
												// surrounded
		this.mcpSpatial = mcpSpatial;
	}

	@ScheduledMethod(start = 2, interval = 1, pick = 1)
	public void ecmRestructure() {
		Context context = ContextUtils.getContext(this);
		fillHoles();
		removeLowCollECM();
		cleanEdge();
		removeExtensions();
	}

	@ScheduledMethod(start = 2, interval = 1)
	public void breakdownECM() {
		// breakdown ECM with MMPs
		Context context = ContextUtils.getContext(this);
		double[] growthFactors = GrowthFactors.getGrowthFactors();
		// if there are MMPs-- chance of breaking down ECM
		if (Fiber.macDepletion == 0 && growthFactors[4] / Fiber.origFiberNumber > 4
				&& this.getCollagen() > 1.2 * Fiber.mdxBaseCollagen
				&& RandomHelper.nextIntFromTo(0, (int) (30 - (growthFactors[4] / Fiber.origFiberNumber))) < 1) {
			// chance of breakdown scales with baseline amount of collagen
			this.setCollagen(getCollagen() - .1);
		}
	}

	public void fillHoles() {
		// If the ECM is surrounded by fibers on 4 sides delete the hole and make the
		// fiber smaller
		Context context = ContextUtils.getContext(this);
		List<Object> ecms = getECM(context);
		for (Object neighbors : ecms) { // for each ecm- get 4 neighbors
			int holeCheck = 0; // reset for each new ECM
			int fiberNumber = -1;
			MooreQuery<Object> query = new MooreQuery(grid, neighbors, 1, 1); //
			Iterable<Object> iter = query.query();
			for (Object neighborCheck : iter) {
				if (neighborCheck instanceof Fiber) {
					fiberNumber = ((Fiber) neighborCheck).getFiberNumber();
					holeCheck = holeCheck + 1; // counter for number of surrounding fibers
				}
			}
			if (holeCheck >= 6) {
				// Then this is ecm surrounded by fibers
				// 1. change to fiber 2. change an exterior edge to be ecm
				MooreQuery<Object> queryOtherFiber = new MooreQuery(grid, neighbors, 1, 1); // get neighbors
				Iterable<Object> iterOtherFiber = queryOtherFiber.query();
				int bordersFiber = 0;
				for (Object otherFibers : iterOtherFiber) {
					if (otherFibers instanceof Fiber && ((Fiber) otherFibers).getFiberNumber() != fiberNumber) {
						bordersFiber = 1; // check if it borders a fiber
					}
				}
				if (bordersFiber == 0) { // the ecm is not bordering fibers with two different fiber numbers
					GridPoint pt = grid.getLocation((ECM) neighbors);
					context.remove((ECM) neighbors);
					Fiber newFiberElem = new Fiber(mcpSpatial, grid, 0, fiberNumber, 0, 0, 0, 0); // change to a fiber
					context.add((Fiber) newFiberElem); // change to a fiber
					grid.moveTo(newFiberElem, pt.getX(), pt.getY());
					// Now have to remove a fiber from the exterior
					List<Object> fiberBorders = ((Fiber) newFiberElem).getFiberBorder(fiberNumber, context); // get the
																												// border
					int index = RandomHelper.nextIntFromTo(0, fiberBorders.size() - 1);
					Object randomBorder = fiberBorders.get(index);
					GridPoint ptRandom = grid.getLocation(randomBorder);
					context.remove(randomBorder); // remove the fiber from the border and make it ECM
					ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to a ecm
					context.add((ECM) newECM);
					grid.moveTo(newECM, ptRandom.getX(), ptRandom.getY());
				}

			}
		}
	}

	public void removeExtensions() {
		// Make sure no fiber pieces are sticking out in long extension
		Context context = ContextUtils.getContext(this);
		List<Object> fibers = Fiber.getFiberElems(context);
		if (fibers != null) {
			for (Object neighbors2 : fibers) { // loop through all fibers, and check if the fiber is a long extension
				int extensionCheck = 0; // reset for each new ECM
				int fiberNumber = ((Fiber) neighbors2).getFiberNumber(); // get the fiber number each time
				MooreQuery<Object> query2 = new MooreQuery(grid, neighbors2, 1, 1); //
				Iterable<Object> iter2 = query2.query();
				for (Object neighborCheck : iter2) {
					if (neighborCheck instanceof ECM) {
						extensionCheck = extensionCheck + 1;
					}
				}
				if (extensionCheck >= 6 && ((Fiber) neighbors2).getElemInFiber(fiberNumber, context).size() > 6) {
					// If the fiber is only 6 elems or smaller- do not get rid of it (keep as a marker in order to regrow)
					// Then this is fiber surrounded by ecm on 7 sides
					// 1. change to ecm 2. change an exterior edge to be fiber
					GridPoint pt2 = grid.getLocation((Fiber) neighbors2);
					// List<Object> fiberBorders2 = ((Fiber) neighbors2).getFiberBorder(fiberNumber,
					// context); // get the border fibers before deleting
					context.remove((Fiber) neighbors2); // remove from context and change to ECM
					ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to ECM
					context.add((ECM) newECM);
					grid.moveTo(newECM, pt2.getX(), pt2.getY());
					// >>>>>>>>Now add a fiber at the edge
					int removeExt = 1; // this marker enforces that the added element will just replace collagen
					double newCollagenVal = ((Fiber) neighbors2).addFiberElem(fiberNumber, context, removeExt);
					// now place the collagen from the removed piece on the added ECM that replaces
					// the removed fiber
					if (newCollagenVal > 0.0) {
						((ECM) newECM).setCollagen(newCollagenVal);
					}
				}
			}
		}
	}

	@ScheduledMethod(start = 2, interval = 24, pick = 1) // only do this every 5 time steps for computational speed
	public void maintainBorderECM() {
		// Keep ECM between the edge and outside fibers-- (no fibers at border of abm)
		Context context = ContextUtils.getContext(this);
		List<Object> ecmEdge = defineEdge(context, grid); // get edge of ECM
		for (Object ecmEdgeTemp : ecmEdge) {
			int fiberNeighborCheck = 0;
			VNQuery<Object> query = new VNQuery(grid, ecmEdgeTemp, 1, 1); // get 4 neighbors
			Iterable<Object> iter = query.query();
			for (Object neighbors : iter) {
				if (neighbors instanceof Fiber) {
					fiberNeighborCheck = 1; // if the ecm edge borders any fibers-- mark with a 1
				}
			}
			if (fiberNeighborCheck == 1) {
				// if it does border a fiber, add an ecm in an open nearby location
				GridPoint openCheck = grid.getLocation(ecmEdgeTemp);
				int[] openCheckX = new int[8];
				int[] openCheckY = new int[8];
				int indexTemp = 0;
				// find open elements in the abm with no agents there-- check from edges of abm
				// inwards
				if (openCheck.getX() - 1 > 0 && openCheck.getY() - 1 > 0 && openCheck.getX() + 1 < muscleBuilder.gridx
						&& openCheck.getY() + 1 < muscleBuilder.gridy) {
					if (!grid.getObjectsAt(openCheck.getX() - 1, openCheck.getY() + 1).iterator().hasNext()) { // if
																												// there
																												// is no
																												// object
																												// -->
																												// open
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() - 1;
						openCheckY[indexTemp - 1] = openCheck.getY() + 1;
					}
					if (!grid.getObjectsAt(openCheck.getX(), openCheck.getY() + 1).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX();
						openCheckY[indexTemp - 1] = openCheck.getY() + 1;
					}
					if (!grid.getObjectsAt(openCheck.getX() + 1, openCheck.getY() + 1).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() + 1;
						openCheckY[indexTemp - 1] = openCheck.getY() + 1;
					}
					if (!grid.getObjectsAt(openCheck.getX() + 1, openCheck.getY()).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() + 1;
						openCheckY[indexTemp - 1] = openCheck.getY();
					}
					if (!grid.getObjectsAt(openCheck.getX() + 1, openCheck.getY() - 1).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() + 1;
						openCheckY[indexTemp - 1] = openCheck.getY() - 1;
					}
					if (!grid.getObjectsAt(openCheck.getX(), openCheck.getY() - 1).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX();
						openCheckY[indexTemp - 1] = openCheck.getY() - 1;
						;
					}
					if (!grid.getObjectsAt(openCheck.getX() - 1, openCheck.getY() - 1).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() - 1;
						openCheckY[indexTemp - 1] = openCheck.getY() - 1;
					}
					if (!grid.getObjectsAt(openCheck.getX() - 1, openCheck.getY()).iterator().hasNext()) {
						indexTemp++;
						openCheckX[indexTemp - 1] = openCheck.getX() - 1;
						openCheckY[indexTemp - 1] = openCheck.getY();
					}
				}
				int randomInt = RandomHelper.nextIntFromTo(0, indexTemp - 1); // indexTemp = number of open spaces
				// Add ECM
				ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to ECM
				context.add((ECM) newECM);
				grid.moveTo(newECM, openCheckX[randomInt], openCheckY[randomInt]);
			}
		}
	}

	public void removeLowCollECM() {
		// Chance of compacting 2 ecm elements with very low collagen
		Context context = ContextUtils.getContext(this);
		List<Object> ecms = ECM.getECM(context);
		// double lowCollTemp = .4;
		double lowCollTemp = 1;
		if (ecms != null) {
			for (Object ecmTemp : ecms) {
				if (ecmTemp instanceof ECM && ((ECM) ecmTemp).getCollagen() < lowCollTemp
						&& RandomHelper.nextIntFromTo(0, 100) < 1) {
					// If the ecm has very low collagen- there is a chance it will be removed and
					// compacted
					GridPoint pt = grid.getLocation((ECM) ecmTemp);
					// Put the collagen on the neighbor
					MooreQuery<Object> query = new MooreQuery(grid, ecmTemp, 1, 1); //
					Iterable<Object> iter = query.query();
					List<Object> ecmNeighbors = new ArrayList();
					for (Object neighbors : iter) {
						if (neighbors instanceof ECM) {
							ecmNeighbors.add(neighbors);
						}
					}
					if (ecmNeighbors.size() > 1) {
						// If there is no neighboring ECM it will stay ECM
						int randomInt = RandomHelper.nextIntFromTo(0, ecmNeighbors.size() - 1);
						Object ecmRandom = ecmNeighbors.get(randomInt);
						// Place current ecm on the neighbor ecm
						((ECM) ecmRandom).setCollagen(((ECM) ecmTemp).getCollagen() + ((ECM) ecmRandom).getCollagen());
						context.remove((ECM) ecmTemp); // remove ECM
						// Need to replace the hole with a neighbor ECM
						restructureECM(pt);
					}
				}
			}
		}

	}

	public void moveByOne(int oldx, int oldy, int newx, int newy) {
		while (grid.getObjectAt(oldx, oldy) != null) { // move all the agents
			grid.moveTo(grid.getObjectAt(oldx, oldy), newx, newy);
		}
	}

	public void restructureECM(GridPoint pt) {
		// Move the entire x or y direction to close the hole
		Context context = ContextUtils.getContext(this);
		// Find what quadrant it is in
		if (pt.getY() > muscleBuilder.gridy / 2.) { // TOP
			if (pt.getX() > muscleBuilder.gridx / 2.) { // TOP RIGHT
				if ((muscleBuilder.gridx - pt.getX()) < (muscleBuilder.gridy - pt.getY())) {
					// Closer in the x direction: move the entire row in the x direction over 1
					for (int i = 1; i < (muscleBuilder.gridx - pt.getX()); i++) {
						int oldx = pt.getX() + i;
						int oldy = pt.getY();
						int newx = pt.getX() + i - 1;
						int newy = oldy;
						moveByOne(oldx, oldy, newx, newy);
					}
				} else {
					// Closer in the y direction; move the entire row in the y direction over 1
					for (int i = 1; i < (muscleBuilder.gridy - pt.getY()); i++) {
						int oldx = pt.getX();
						int oldy = pt.getY() + i;
						int newx = oldx;
						int newy = pt.getY() + i - 1;
						moveByOne(oldx, oldy, newx, newy);
					}
				}
			} else { // TOP LEFT
				if (pt.getX() < (muscleBuilder.gridy - pt.getY())) {
					// Closer in the x direction: move the entire row in the x direction over 1
					for (int i = 1; i < pt.getX(); i++) {
						int oldx = pt.getX() - i;
						int oldy = pt.getY();
						int newx = pt.getX() - i + 1;
						int newy = oldy;
						moveByOne(oldx, oldy, newx, newy);
					}
				} else {
					// Closer in the y direction; move the entire row in the y direction over 1
					for (int i = 1; i < (muscleBuilder.gridy - pt.getY()); i++) {
						int oldx = pt.getX();
						int oldy = pt.getY() + i;
						int newx = oldx;
						int newy = pt.getY() + i - 1;
						moveByOne(oldx, oldy, newx, newy);
					}
				}
			}
		} else { // BOTTOM
			if (pt.getX() > muscleBuilder.gridx / 2.) { // BOTTOM RIGHT
				if ((muscleBuilder.gridx - pt.getX()) < pt.getY()) {
					// Closer in the x direction: move the entire row in the x direction over 1
					for (int i = 1; i < (muscleBuilder.gridx - pt.getX()); i++) {
						int oldx = pt.getX() + i;
						int oldy = pt.getY();
						int newx = pt.getX() + i - 1;
						int newy = oldy;
						moveByOne(oldx, oldy, newx, newy);
					}
				} else {
					// Closer in the y direction; move the entire row in the y direction over 1
					for (int i = 1; i < pt.getY() - 1; i++) {
						int oldx = pt.getX();
						int oldy = pt.getY() - i;
						int newx = oldx;
						int newy = pt.getY() - i + 1;
						moveByOne(oldx, oldy, newx, newy);
					}
				}
			} else { // BOTTOM LEFT
				if (pt.getX() < pt.getY()) {
					// Closer in the x direction: move the entire row in the x direction over 1
					for (int i = 1; i < pt.getX() - 1; i++) {
						int oldx = pt.getX() - i;
						int oldy = pt.getY();
						int newx = pt.getX() - i + 1;
						int newy = oldy;
						moveByOne(oldx, oldy, newx, newy);
					}
				} else {
					// Closer in the y direction; move the entire row in the y direction over 1
					for (int i = 1; i < pt.getY() - 1; i++) {
						int oldx = pt.getX();
						int oldy = pt.getY() - i;
						int newx = oldx;
						int newy = pt.getY() - i + 1;
						moveByOne(oldx, oldy, newx, newy);
					}
				}
			}
		}
	}

	public static List<Object> defineEdge(Context<Object> context, Grid<Object> grid) {
		// Get all the ECM
		List<Object> ecms = getECM(context); // create a list of all the ecm agents
		List<Object> ecmEdge = new ArrayList<Object>(); // get the ecm edge
		// Check neighbors
		if (ecms != null) {
			for (Object ecmTemp : ecms) {
				// get all the ecm neighbors
				MooreQuery<Object> query = new MooreQuery(grid, ecmTemp, 1, 1);
				Iterable<Object> iter = query.query();
				List<Object> neighborsList = new ArrayList();
				for (Object neighborTemp : iter) {
					if (neighborTemp instanceof Necrosis || neighborTemp instanceof ECM
							|| neighborTemp instanceof Fiber) {
						neighborsList.add(neighborTemp); // determine how many neighbors there are
					}
				}
				if (neighborsList.size() < 8) {
					((ECM) ecmTemp).setEcmEdge(1);
					ecmEdge.add((ECM) ecmTemp); // add to the ecm edge list
				} else {
					((ECM) ecmTemp).setEcmEdge(0);
				}
			}
		}
		return ecmEdge;
	}

	public List<Object> cleanEdge() {
		// Clean up the ecm edges after the restructuring so there are not long
		// protrusions/holes
		// Get the edge of the grid
		Context context = ContextUtils.getContext(this);
		List<Object> ecms = new ArrayList<Object>(); // create a list of all the ecm agents
		List<Object> ecmEdge = new ArrayList<Object>();
		if (context != null) {
			for (Object obj : context) {
				if (obj instanceof ECM) {
					ecms.add(obj);
				}
			}
		}
		for (Object ecmTemp : ecms) {
			// get all the ecm neighbors
			MooreQuery<Object> query = new MooreQuery(grid, ecmTemp, 1, 1); //
			Iterable<Object> iter = query.query();
			List<Object> neighborsList = new ArrayList();
			for (Object neighborTemp : iter) {
				if (neighborTemp instanceof Necrosis || neighborTemp instanceof ECM || neighborTemp instanceof Fiber) {
					neighborsList.add(neighborTemp); // determine how many neighbors there are
				}
			}
			if (neighborsList.size() < 8) {
				((ECM) ecmTemp).setEcmEdge(1);
				ecmEdge.add((ECM) ecmTemp);
				// Also check if it is a protrusion
				if (neighborsList.size() < 3) {
					// remove and restructure the protrusion
					// a fiber edge that has no ecm outside of it (at edge)
					context.remove((ECM) ecmTemp);
					fillEdgeHole();
				}
			} else {
				((ECM) ecmTemp).setEcmEdge(0);
			}
		}
		return ecmEdge;
	}

	public void fillEdgeHole() {
		// Fill holes in edge
		Context context = ContextUtils.getContext(this);
		List<Object> fibers = Fiber.getFiberElems(context);
		List<Object> fiberEdgeNoECM = new ArrayList();
		if (fibers != null) {
			for (Object fiberTemp : fibers) {
				// get all the ecm neighbors
				MooreQuery<Object> query = new MooreQuery(grid, fiberTemp, 1, 1); //
				Iterable<Object> iter = query.query();
				List<Object> neighborsList = new ArrayList();
				for (Object neighborTemp : iter) {
					neighborsList.add(neighborTemp); // determine how many neighbors there are
				}
				if (neighborsList.size() < 8) {
					fiberEdgeNoECM.add(fiberTemp);
				}
			}
		}
		if (fiberEdgeNoECM.size() > 0) {
			// Choose a random fiber
			int random = RandomHelper.nextIntFromTo(0, fiberEdgeNoECM.size() - 1);
			Object fibRandom = fiberEdgeNoECM.get(random);
			// Put the ecm in the open space surrounding the fiber
			// Find empty site:
			GridPoint pt = grid.getLocation(fibRandom);
			int[] emptySites = findEmptySites(pt, grid);
			if (emptySites[0] != -1) {
				// Then move the collagen to this empty site
				ECM newECM = new ECM(mcpSpatial, grid, 1, 1, 0); // change to ECM
				context.add((ECM) newECM);
				grid.moveTo(newECM, emptySites[0], emptySites[1]);
			}
		} else {
			// if there are no fiber edges that need to be refilled-- move to ECM edge
			// neighbor
			List<Object> ecmEdgeElems = getECMEdgeElems(context);
			if (ecmEdgeElems != null) {
				int random2 = RandomHelper.nextIntFromTo(0, ecmEdgeElems.size() - 1);
				Object ecmRandom = ecmEdgeElems.get(random2);
				GridPoint pt = grid.getLocation(ecmRandom);
				int[] emptySites = findEmptySites(pt, grid);
				if (emptySites[0] != -1) {
					// Then move the collagen to this empty site
					ECM newECM = new ECM(mcpSpatial, grid, 1, 0, 0); // change to ECM
					context.add((ECM) newECM);
					grid.moveTo(newECM, emptySites[0], emptySites[1]);
				}
			}
		}
	}

	public void ecmBehaviors() {
		// ECM between fibers should be maintained
	}

	public static List<Object> getECMEdgeElems(Context<Object> context) {
		// get ecm elems (agents) in the edge
		List<Object> ecms = ECM.getECM(context);
		List<Object> ecmEdgeElem = new ArrayList();
		if (ecms != null) {
			for (Object ecmTemp : ecms) {
				if (((ECM) ecmTemp).ecmEdge == 1) {
					ecmEdgeElem.add((ECM) ecmTemp);
				}
			}
			return ecmEdgeElem;
		} else
			return null;
	}

	public static int[] findEmptySites(GridPoint pt, Grid<Object> grid) {
		int[] emptySites = new int[2];
		int newPtx = -1;
		int newPty = -1;
		if (pt.getX() - 1 > 0 && pt.getY() - 1 > 0 && pt.getX() + 1 < muscleBuilder.gridx
				&& pt.getY() + 1 < muscleBuilder.gridy) {
			if (!grid.getObjectsAt(pt.getX() - 1, pt.getY() + 1).iterator().hasNext()) {
				newPtx = pt.getX() - 1;
				newPty = pt.getY() + 1;
			} else if (!grid.getObjectsAt(pt.getX(), pt.getY() + 1).iterator().hasNext()) {
				newPtx = pt.getX();
				newPty = pt.getY() + 1;
			} else if (!grid.getObjectsAt(pt.getX() + 1, pt.getY() + 1).iterator().hasNext()) {
				newPtx = pt.getX() + 1;
				newPty = pt.getY() + 1;
			} else if (!grid.getObjectsAt(pt.getX() + 1, pt.getY()).iterator().hasNext()) {
				newPtx = pt.getX() + 1;
				newPty = pt.getY();
			} else if (!grid.getObjectsAt(pt.getX() + 1, pt.getY() - 1).iterator().hasNext()) {
				newPtx = pt.getX() + 1;
				newPty = pt.getY() - 1;
			} else if (!grid.getObjectsAt(pt.getX(), pt.getY() - 1).iterator().hasNext()) {
				newPtx = pt.getX();
				newPty = pt.getY() - 1;
			} else if (!grid.getObjectsAt(pt.getX() - 1, pt.getY() - 1).iterator().hasNext()) {
				newPtx = pt.getX() - 1;
				newPty = pt.getY() - 1;
			} else if (!grid.getObjectsAt(pt.getX() - 1, pt.getY()).iterator().hasNext()) {
				newPtx = pt.getX() - 1;
				newPty = pt.getY();
			}
		}
		emptySites[0] = newPtx;
		emptySites[1] = newPty;
		return emptySites;
	}

	public static List<Object> getECM(Context<Object> context) {
		// Get a list of all the ecm agents to place cells on the ECM
		List<Object> ecms = new ArrayList<Object>(); // create a list of all the ecm agents
		if (context != null) {
			for (Object obj : context) {
				if (obj instanceof ECM) {
					ecms.add(obj);
				}
			}
			return ecms;
		} else
			return null;
	}

	public static double getTotalCollagenAmt(Context<Object> context) {
		// add the total amount of collagen in all the ecm agents/elements
		double totalCollagen = 0;
		List<Object> ecms = getECM(context);
		for (Object obj : ecms) {
			totalCollagen = totalCollagen + ((ECM) obj).getCollagen();
		}
		return totalCollagen;
	}

	public int getEcmEdge() {
		return this.ecmEdge;
	}

	public void setEcmEdge(int ecmEdge) {
		this.ecmEdge = ecmEdge;
	}

	public double getCollagen() {
		return collagen;
	}

	public void setCollagen(double collagen) {
		this.collagen = collagen;
	}

	public int getSameFiberBorder() {
		return this.sameFiberBorder;
	}

	public void setSameFiberBorder(int sameFiberBorder) {
		this.sameFiberBorder = sameFiberBorder;
	}

}
