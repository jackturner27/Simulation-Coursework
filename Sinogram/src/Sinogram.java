
import java.util.Arrays;

import java.awt.*;
import javax.swing.*;

public class Sinogram {

	static final int N = 512;

	static final int CELL_SIZE = 1;

	static final double SCALE = 0.0045; // think of a better way to
										// parametrize this later...

	static final int CUTOFF = N / 4; // in ramp filter

	static final float GREY_SCALE_LO = 0.95f, GREY_SCALE_HI = 1.05f;
	// Clipping, for display only. See for example Figure 1 in:
	// http://bigwww.epfl.ch/thevenaz/shepplogan/

	public static void main(String[] args) {

		double[][] density = new double[N][N];

		for (int i = 0; i < N; i++) {
			double x = SCALE * (i - N / 2);
			for (int j = 0; j < N; j++) {
				double y = SCALE * (j - N / 2);

				density[i][j] = sheppLoganPhantom(x, y);
			}
		}

		DisplayDensity display1 = new DisplayDensity(density, N, "Source Model", GREY_SCALE_LO, GREY_SCALE_HI);

		// Radon tranform of density (as measured by detectors):

		double[][] sinogram = new double[N][N];

		for (int iTheta = 0; iTheta < N; iTheta++) {
			double theta = (Math.PI * iTheta) / N;
			double cos = Math.cos(theta);
			double sin = Math.sin(theta);
			for (int iR = 0; iR < N; iR++) {
				double r = SCALE * (iR - N / 2);
				double sum = 0;
				for (int iS = 0; iS < N; iS++) {
					double s = SCALE * (iS - N / 2);
					double x = r * cos + s * sin;
					double y = r * sin - s * cos;
					sum += sheppLoganPhantom(x, y);
				}
				sinogram[iTheta][iR] = sum;
			}
		}

		DisplayDensity display2 = new DisplayDensity(sinogram, N, "Sinogram");

		// inferred integral of density points (actually sum of density
		// points, here) for laternormalization of reconstruction

		double normDensity = norm1(sinogram[0]);

		// ... Insert sinogram filtering code here! ...

		double[][] sinogramFTRe = new double[N][N], sinogramFTIm = new double[N][N];
		for (int iTheta = 0; iTheta < N; iTheta++) {
			for (int iR = 0; iR < N; iR++) {
				sinogramFTRe[iTheta][iR] = sinogram[iTheta][iR];
			}
		}

		for (int iTheta = 0; iTheta < N; iTheta++) {
			// ... do 1D FFT on a row ...
			FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], 1);
		}

		DisplaySinogramFT display3 = new DisplaySinogramFT(sinogramFTRe, sinogramFTIm, N,
				"Sinogram radial Fourier Transform");

		for (int iTheta = 0; iTheta < N; iTheta++) {
			for (int iK = 0; iK < N; iK++) {
				int kSigned = iK <= N / 2 ? iK : iK - N;
				// ... multiply Sinogram FT by abs(kSigned) ...
				
				double arg = Math.abs(kSigned) * Math.cos(Math.PI * kSigned / ( 2 * CUTOFF ));
				
				sinogramFTRe[iTheta][iK] *= Math.abs(kSigned);
				sinogramFTIm[iTheta][iK] *= Math.abs(kSigned);
			}
		}

		for (int iTheta = 0; iTheta < N; iTheta++) {
			// ... do inverse 1D FFT on a row ...
			FFT.fft1d(sinogramFTRe[iTheta], sinogramFTIm[iTheta], -1);
		}

		DisplayDensity display4 = new DisplayDensity(sinogramFTRe, N, "Filtered sinogram");

		double[][] backProjection = new double[N][N];
		backProject(backProjection, sinogramFTRe);

		// Normalize reconstruction, to have same sum as inferred for
		// original density

		double factor = normDensity / norm2(backProjection);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				backProjection[i][j] *= factor;
			}
		}

		DisplayDensity display5 = new DisplayDensity(backProjection, N, "Back projected sinogram",
                GREY_SCALE_LO, GREY_SCALE_HI);
	}

	static void backProject(double[][] projection, double[][] sinogram) {

		// Back Projection operation on sinogram

		for (int i = 0; i < N; i++) {
			double x = SCALE * (i - N / 2);
			for (int j = 0; j < N; j++) {
				double y = SCALE * (j - N / 2);

				double sum = 0;
				for (int iTheta = 0; iTheta < N; iTheta++) {
					double theta = (Math.PI * iTheta) / N;
					double cos = Math.cos(theta);
					double sin = Math.sin(theta);

					double r = x * cos + y * sin;

					double rBox = N / 2 + r / SCALE;

					if (rBox < 0)
						continue; // assume centred object, with
									// support radius < N/2

					int iR = (int) rBox;

					double offR = rBox - iR;
					int iPlusR = iR + 1;

					if (iPlusR >= N)
						continue; // ditto.

					// linear interpolation
					double sinogramVal = (1 - offR) * sinogram[iTheta][iR] + offR * sinogram[iTheta][iPlusR];
					sum += sinogramVal;
				}
				projection[i][j] = sum;
			}
		}
	}

	// Shepp-Logan Phantom:
	//
	// https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom

	static final Ellipse[] sheppLoganEllipses = { new Ellipse(0.0, 0.0, 0.69, 0.92, 0, 2.0),
			new Ellipse(0.0, -0.0184, 0.6624, 0.874, 0, -0.98), new Ellipse(0.22, 0, 0.11, 0.31, -18.0, -0.02),
			new Ellipse(-0.22, 0, 0.16, 0.41, 18.0, -0.02), new Ellipse(0, 0.35, 0.21, 0.25, 0, 0.01),
			new Ellipse(0, 0.1, 0.046, 0.046, 0, 0.01), new Ellipse(0, -0.1, 0.046, 0.046, 0, 0.01),
			new Ellipse(-0.08, -0.605, 0.046, 0.023, 0, 0.01), new Ellipse(0, -0.605, 0.023, 0.023, 0, 0.01),
			new Ellipse(0.06, -0.605, 0.023, 0.046, 0, 0.01), };

	static double sheppLoganPhantom(double x, double y) {

		double total = 0;
		for (Ellipse ellipse : sheppLoganEllipses) {
			total += ellipse.localDensity(x, y);
		}
		return total;
	}

	static class Ellipse {

		double centreX;
		double centreY;
		double major;
		double minor;
		double theta;
		double density;
		double cos, sin;

		Ellipse(double centreX, double centreY, double major, double minor, double theta, double density) {

			this.centreX = centreX;
			this.centreY = centreY;
			this.major = major;
			this.minor = minor;
			this.theta = theta;
			if (theta == 0) {
				cos = 1;
				sin = 0;
			} else {
				double rad = Math.PI * theta / 180;
				cos = Math.cos(rad);
				sin = Math.sin(rad);
			}
			this.density = density;
		}

		double localDensity(double x, double y) {

			double xOff, yOff;
			xOff = x - centreX;
			yOff = y - centreY;

			double xRot, yRot;
			if (theta == 0) {
				xRot = xOff;
				yRot = yOff;
			} else {
				// Rotate so x/y aligned with major/minor axes.
				xRot = cos * xOff - sin * yOff;
				yRot = sin * xOff + cos * yOff;
			}
			double xNorm = xRot / major;
			double yNorm = yRot / minor;
			if (xNorm * xNorm + yNorm * yNorm < 1) {
				return density;
			} else {
				return 0;
			}
		}
	}

	static double norm1(double[] density) {

		double norm = 0;
		for (int i = 0; i < N; i++) {
			norm += density[i];
		}
		return norm;
	}

	static double norm2(double[][] density) {

		double norm = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (density[i][j] > 0) {
					norm += density[i][j];
				}
			}
		}
		return norm;
	}

}
