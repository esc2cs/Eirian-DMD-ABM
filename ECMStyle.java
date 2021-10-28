/**
 * 
 */
package muscle;

import java.awt.Color;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.ShapeFactory2D;
import saf.v3d.scene.VSpatial;

/**
 * @author Kelley Virgilio
 *
 */
public class ECMStyle extends DefaultStyleOGL2D {

	private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory) {
		this.shapeFactory = factory;
	}

	@Override
	public Color getColor(Object o) {

		// return new Color(196,196,196); // ECM all one color
		return new Color(207, 174, 242); // ECM all one color
	}

	@Override
	public VSpatial getVSpatial(Object agent, VSpatial spatial) {
		if (spatial == null) {
			spatial = shapeFactory.createRectangle(15, 15);
		}
		return spatial;
	}
}
