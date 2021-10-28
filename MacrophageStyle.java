package muscle;

import java.awt.Color;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.ShapeFactory2D;
import saf.v3d.scene.VSpatial;

/**
 * @author Kelley Virgilio
 *
 */
 
public class MacrophageStyle extends DefaultStyleOGL2D{

private ShapeFactory2D shapeFactory;

	@Override
	public void init(ShapeFactory2D factory){
		this.shapeFactory = factory;
	}
	
	@Override
	public Color getColor(Object obj) {
		if(((Macrophage) obj).getPhenotype() == 1) {
			return new Color(230, 121, 132);//a light coral (orangy pink) color
		}
		else if(((Macrophage) obj).getPhenotype() == 2) {
			return new Color(214, 208, 222);//light gray
		}
		else {
			return new Color(170, 170, 170);//dark gray
		}
	}
	
	@Override
	  public VSpatial getVSpatial(Object agent, VSpatial spatial) {
	    if (spatial == null) {
	      //spatial = shapeFactory.createRectangle(25,25);
	      spatial = shapeFactory.createRectangle(30,30);
	    }
	    return spatial;
	  }

}