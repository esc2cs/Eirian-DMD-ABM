package muscle;

import java.awt.Color;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.ShapeFactory2D;
import saf.v3d.scene.VSpatial;

/**
 * @author Kelley Virgilio
 *
 */

public class FibroblastStyle extends DefaultStyleOGL2D {

	private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory) {
		this.shapeFactory = factory;
	}

	@Override
	public Color getColor(Object o) {
		if (((Fibroblast) o).getResident() == 1) {
			return new Color(252, 230, 167); // not active- resident- light yellow
		} else if (((Fibroblast) o).getResident() == 0
				&& ((Fibroblast) o).getPhenotype() >= Fibroblast.myofibroblastSwitch) {
			return new Color(140, 63, 7); // myofibroblast- dark orange
		}
		else if (((Fibroblast) o).getResident() == 0
				&& ((Fibroblast) o).getPhenotype() <= Fibroblast.fapSwitch) {
			return new Color(156, 0, 0); // FAP- dark red
		}

		else {
			return new Color(237, 142, 74); // regular active fibroblast- light orange
		}

	}

	@Override
	public VSpatial getVSpatial(Object agent, VSpatial spatial) {
		if (spatial == null) {
			// spatial = shapeFactory.createRectangle(25,25);
			spatial = shapeFactory.createRectangle(30, 30);
		}
		return spatial;
	}

}
