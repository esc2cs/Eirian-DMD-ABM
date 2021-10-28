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
public class SSCStyle extends DefaultStyleOGL2D {

	private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory) {
		this.shapeFactory = factory;
	}

	@Override
	public Color getColor(Object o) {

		// ACTIVATED SSC- dark purple
		if (((SSC) o).getActive() >= 4 && ((SSC) o).getDiff() == 0 && ((SSC) o).getCommitted() != 2) { // if active and
																										// NOT
																										// differentiated,
																										// NOT committed
																										// = activated
																										// SSC
			return new Color(157, 70, 219); // darker purple
		}
		// MYOBLAST - blue
		else if (((SSC) o).getActive() >= 4 && ((SSC) o).getDiff() == 0 && ((SSC) o).getCommitted() == 2) { // myoblast
			return new Color(0, 102, 204); // blue
		}
		// MYOCYTE- navy blue
		else if (((SSC) o).getActive() >= 4 && ((SSC) o).getDiff() == 2) {
			return new Color(39, 30, 201); // myocyte that is differentiating or not on a fiber- dark blue/navy
		}
		// MYOCYTE ON EDGE --> REGROWING MUSCLE- hot pink
		else if (((SSC) o).getActive() >= 4 && ((SSC) o).getDiff() == 3) {
			return new Color(255, 0, 127); // fused differentiated myocyte- hot pink
		}
		// QUIESCENT SSC OR MYOBLAST- light green
		else
			return new Color(131, 232, 128); // light green when quiescent
	}

	@Override
	public VSpatial getVSpatial(Object agent, VSpatial spatial) {
		if (spatial == null) {
			// spatial = shapeFactory.createCircle(15,15);
			spatial = shapeFactory.createCircle(20, 20);
		}
		return spatial;
	}
}
