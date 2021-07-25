public class FFTImageFiltering {

    public static int N = 256 ;

    public static void main(String [] args) throws Exception {

        double [] [] X = new double [N] [N] ;
        ReadPGM.read(X, "wolf.pgm", N) ;

        DisplayDensity display =
                new DisplayDensity(X, N, "Original Image") ;

        // create array for in-place FFT, and copy original data to it
        double [] [] CRe = new double [N] [N], CIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                CRe [k] [l] = X [k] [l] ;
            }
        }

        fft2d(CRe, CIm, 1) ;  // Fourier transform

        Display2dFT display2 =
                new Display2dFT(CRe, CIm, N, "Discrete FT") ;

        // create array for in-place inverse FFT, and copy FT to it
        double [] [] reconRe = new double [N] [N],
                     reconIm = new double [N] [N] ;
        for(int k = 0 ; k < N ; k++) {
            for(int l = 0 ; l < N ; l++) {
                reconRe [k] [l] = CRe [k] [l] ;
                reconIm [k] [l] = CIm [k] [l] ;
            }
        }

        fft2d(reconRe, reconIm, -1) ;  // Inverse Fourier transform

        DisplayDensity display3 =
                new DisplayDensity(reconRe, N, "Reconstructed Image") ;
    }

    
    static void transpose(double [] [] a) {

    	double [][] b = new double[N][N];
    	
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < i ; j++) {
                 //... swap values in a [i] [j] and a [j] [i] elements ...
            	b[i][j] = a[j][i];
            }
        }
    }
    
    static void fft2d(double [] [] re, double [] [] im, int isgn) {

        // For simplicity, assume square arrays

        //... fft1d on all rows of re, im ...
    	for(int i = 0; i < N; i++) {
    		FFT.fft1d(re[i], im[i], isgn);
    	}
    	
        transpose(re) ;
        transpose(im) ;

        //... fft1d on all rows of re, im ...
        for(int i = 0; i < N; i++) {
    		FFT.fft1d(re[i], im[i], isgn);
    	}

        transpose(re) ;
        transpose(im) ;
    }

}