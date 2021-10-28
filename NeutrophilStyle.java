/**
 * 
 */
package muscle;

import java.awt.Color;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.ShapeFactory2D;
import saf.v3d.scene.VSpatial;

/**
 * @author kbc7vn
 *
 */
public class NeutrophilStyle extends DefaultStyleOGL2D {

	private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory) {
		this.shapeFactory = factory;
	}

	@Override
	public Color getColor(Object o) {
		if (((Neutrophil) o).getApoptosed() == 0) {
			return new Color(252, 230, 167); // neutrophil is not apoptosed anything- light yellow
		} else {
			return new Color(140, 63, 7); // apoptosed neutrophil - brown
		}

	}

	@Override
	public VSpatial getVSpatial(Object agent, VSpatial spatial) {
		if (spatial == null) {
			// spatial = shapeFactory.createRectangle(25,25);
			spatial = shapeFactory.createCircle(30, 30);
		}
		return spatial;
	}

}
