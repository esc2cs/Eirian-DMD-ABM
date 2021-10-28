package muscle;

import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.valueLayer.BufferedGridValueLayer;
import repast.simphony.valueLayer.GridValueLayer;
import repast.simphony.valueLayer.IGridValueLayer;
import repast.simphony.valueLayer.ValueLayer;
import repast.simphony.valueLayer.ValueLayerDiffuser;

public class MCPContext extends DefaultContext<muscleBuilder> {	

	private ValueLayerDiffuser diffuser;

	/*
	 * (non-Javadoc)
	 * 
	 * @see repast.simphony.context.AbstractContext#addValueLayer(repast.simphony.
	 * valueLayer.ValueLayer)
	 */
	@Override
	public void addValueLayer(ValueLayer valueLayer) {
		// TODO Auto-generated method stub
		super.addValueLayer(valueLayer);
		diffuser = new ValueLayerDiffuser((IGridValueLayer) valueLayer, Math.pow(.5, 1.0 / 5.), 1.0, false); // set evaporation constant to same decay as in GrowthFactors.java
		diffuser.setMinValue(0);
		
	}

	/**
	 * Swaps the buffered heat layers.
	 */
	// priority = -1 so that the heatbugs action occurs first
	@ScheduledMethod(start = 1, interval = 1, priority = -1)
	public void MCPdiffuse() {
	//public void swap() {
		BufferedGridValueLayer grid = (BufferedGridValueLayer) getValueLayer("MCP Layer");
		//GridValueLayer grid = (GridValueLayer) getValueLayer("MCP Layer");
		grid.swap();
		diffuser.diffuse();
	}
	


}
