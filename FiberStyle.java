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
public class FiberStyle extends DefaultStyleOGL2D {

	private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory) {
		this.shapeFactory = factory;
	}

	@Override
	public Color getColor(Object o) {

		if (Fiber.tick > 2) {
			if (((Fiber) o).getBorder() == 1 && ((Fiber) o).getNeedsRepair() == 0) { // if it is boarder fiber change
																						// the color to dark red
				return new Color(168, 0, 0); // dark red
			} else if (((Fiber) o).getBorder() == 1 && ((Fiber) o).getNeedsRepair() == 1) { // if it needs repair
				return new Color(13, 26, 145); // dark blue
			} else
				return new Color(242, 146, 146);
		}

		else
			return new Color(242, 146, 146); // salmon

	}

	@Override
	public VSpatial getVSpatial(Object agent, VSpatial spatial) {
		if (spatial == null) {
			spatial = shapeFactory.createRectangle(15, 15);
		}
		return spatial;
	}
}
