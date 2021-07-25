
//import java.util.Arrays ;

import java.awt.* ;
import javax.swing.* ;

import java.awt.image.BufferedImage ;

import javax.imageio.ImageIO;

import java.io.File ;


public class SheppLoganImage {

    static final int N = 2048 ;

    static final int CELL_SIZE = 1 ;

    //static final double SCALE = 0.018 ;  // think of a better way to
    static final double SCALE = 0.001125 ;  // think of a better way to
                                          // parametrize this later...

    static final float GREY_SCALE_LO = 0.95f, GREY_SCALE_HI = 1.05f ;
        // Clipping, for display only.  See for example Figure 1 in:
        //    http://bigwww.epfl.ch/thevenaz/shepplogan/

    public static void main(String [] args) throws Exception {

        double [] [] density = new double [N] [N] ;

        for(int i = 0 ; i < N ; i++) {
            double x = SCALE * (i - N/2) ;
            for(int j = 0 ; j < N ; j++) {
                double y = SCALE * (j - N/2) ;

                density [i] [j] = sheppLoganPhantom(x, y) ;
            }
        } 

        DisplayDensity display1 =
                new DisplayDensity(density, "Source Model",
                                   GREY_SCALE_LO, GREY_SCALE_HI) ;
                //new DisplayDensity(density, "Source Model") ;

        BufferedImage img = new BufferedImage(N, N,
                                              BufferedImage.TYPE_INT_ARGB) ;

        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                int grey = (int) (255 * density [i] [j] / 2.0 + 0.5) ;
                //int rgb = grey + (grey << 8) + (grey << 16) ;
                Color c = new Color(grey, grey, grey) ;
                img.setRGB(i, j, c.getRGB()) ;
                System.out.println("i = " + i + ", j = " + j + ": " + grey) ;
            }
        } 

        ImageIO.write(img, "PNG", new File("shepp-logan-xlarge.png"));
    }


    // Shepp-Logan Phantom:
    //
    //   https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom

    static final Ellipse [] sheppLoganEllipses = {
        new Ellipse(0.0, 0.0, 0.69, 0.92, 0, 2.0),
        new Ellipse(0.0, -0.0184, 0.6624, 0.874, 0, -0.98),
        new Ellipse(0.22, 0, 0.11, 0.31, -18.0, -0.02),
        new Ellipse(-0.22, 0, 0.16, 0.41, 18.0, -0.02),
        new Ellipse(0, 0.35, 0.21, 0.25, 0, 0.01),
        new Ellipse(0, 0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(0, -0.1, 0.046, 0.046, 0, 0.01),
        new Ellipse(-0.08, -0.605, 0.046, 0.023, 0, 0.01),
        new Ellipse(0, -0.605, 0.023, 0.023, 0, 0.01),
        new Ellipse(0.06, -0.605, 0.023, 0.046, 0, 0.01),
    } ;

    static double sheppLoganPhantom (double x, double y) {

        double total = 0 ;
        for(Ellipse ellipse : sheppLoganEllipses) {
            total += ellipse.localDensity(x, y) ;
        }
        return total ;
    }

    static class Ellipse {

        double centreX ;
        double centreY ;
        double major ;
        double minor ;
        double theta ;
        double density ;
        double cos, sin ;

        Ellipse(double centreX, double centreY,
                double major, double minor, double theta, double density) {

            this.centreX = centreX ;
            this.centreY = centreY ;
            this.major = major ;
            this.minor = minor ;
            this.theta = theta ;
            if(theta == 0) {
                cos = 1 ;
                sin = 0 ;
            }
            else {
                double rad = Math.PI * theta / 180 ;
                cos = Math.cos(rad) ; 
                sin = Math.sin(rad) ; 
            }
            this.density = density;
        }

        double localDensity(double x, double y) {

            double xOff, yOff ;
            xOff = x - centreX ;
            yOff = y - centreY ;

            double xRot, yRot ;
            if(theta == 0) {
                xRot = xOff ;
                yRot = yOff ;
            }
            else {
                // Rotate so x/y aligned with major/minor axes.
                xRot = cos * xOff - sin * yOff ;
                yRot = sin * xOff + cos * yOff ;
            }
            double xNorm = xRot / major ;
            double yNorm = yRot / minor ;
            if(xNorm * xNorm + yNorm * yNorm < 1) {
                return density ;
            }
            else {
                return 0 ;
            }
        }
    }

    static class DisplayDensity extends JPanel {

        //final static int WINDOW_SIZE = N * CELL_SIZE ;

        double [] [] density ;

        double greyScaleLo, greyScaleHi ;

        boolean doScale ;

        DisplayDensity(double [] [] density, String title) {
            this(density, title, Double.MIN_VALUE, Double.MAX_VALUE) ;
        }

        DisplayDensity(double [] [] density, String title,
                       double greyScaleLo, double greyScaleHi) {

            setPreferredSize(new Dimension(CELL_SIZE * N, CELL_SIZE * N)) ;

            JFrame frame = new JFrame(title);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);

            this.density = density ;

            this.greyScaleLo = greyScaleLo ;
            this.greyScaleHi = greyScaleHi ;

            repaint() ;
        }

        public void paintComponent(Graphics g) {

            // If eithe hi or lo clipping unspecified, find min and max
            // data values ;

            double minVal = Double.MAX_VALUE ;
            double maxVal = Double.MIN_VALUE ;

            if(greyScaleLo == Double.MIN_VALUE ||
                greyScaleHi == Double.MAX_VALUE) {
                for(int i = 0 ; i < N ; i++) {
                    for(int j = 0 ; j < N ; j++) {
                        double densityVal = density [i] [j] ;
                        if(densityVal > maxVal) {
                            maxVal = densityVal ;
                        }
                        if(densityVal < minVal) {
                            minVal = densityVal ;
                        }
                    }
                }
            }
            if(greyScaleLo == Double.MIN_VALUE) {
                greyScaleLo = minVal ;
            }
            if(greyScaleHi == Double.MAX_VALUE) {
                greyScaleHi = maxVal ;
            }

            for(int i = 0 ; i < N ; i++) {
                for(int j = 0 ; j < N ; j++) {
                    double intensity = density [i] [N - j - 1] ;

                    float grey ;
                    if(intensity <= greyScaleLo) {
                        grey = 0 ;
                    } else if(intensity >= greyScaleHi) {
                        grey = 1 ;
                    }
                    else {
                        grey = (float) ((intensity - greyScaleLo) /
                                        (greyScaleHi - greyScaleLo)) ;
                    }
                    Color c = new Color(grey, grey, grey) ;
                    g.setColor(c) ;
                    g.fillRect(CELL_SIZE * i, CELL_SIZE * j,
                               CELL_SIZE, CELL_SIZE) ;
                }
            }
        }
    }
}
