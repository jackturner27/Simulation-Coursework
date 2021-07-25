public class FFTImageFiltering {

    public static int N = 256 ;

    public static void main(String [] args) throws Exception {
    	
    	double [] [] X = new double [N] [N] ;
        ReadPGM.read(X, "wolf.pgm", N) ;

        DisplayDensity display =
                new DisplayDensity(X, N, "Original Image") ;

        double starttime = System.currentTimeMillis();
        
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

        
        /*int cutoff = N/8 ;  // for example
        for(int k = 0 ; k < N ; k++) {
            int kSigned = k <= N/2 ? k : k - N ;
            for(int l = 0 ; l < N ; l++) {
                int lSigned = l <= N/2 ? l : l - N ;
                if(Math.abs(kSigned) > cutoff || Math.abs(lSigned) > cutoff) {
                    CRe [k] [l] = 0 ;
                    CIm [k] [l] = 0 ;
                }
            }
        }

        Display2dFT truncated =
                new Display2dFT(CRe, CIm, N, "Truncated FT") ;*/
        
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
        
        double endtime = System.currentTimeMillis();
        System.out.println(endtime - starttime + "ms");
    }
    
    static double [] []  transpose(double [] [] a) {
    	
    	double [] [] b = new double [N] [N]; //create temporary matrix
    	
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < N ; j++) {
                 //... swap values in a [i] [j] and a [j] [i] elements ...
            	b[i][j] = a[j][i];
            }
        }
        return b;
    }
    
    static void fft2d(double [] [] re, double [] [] im, int isgn) {

        // For simplicity, assume square arrays

        //... fft1d on all rows of re, im ...
    	    	
    	for(int i = 0; i < N; i++) {
    		FFT.fft1d(re[i], im[i], isgn);
    	}

        re = transpose(re) ;
        im = transpose(im) ;

        //... fft1d on all rows of re, im ...
        
        for(int i = 0; i < N; i++) {
    		FFT.fft1d(re[i], im[i], isgn);
    	}

        re = transpose(re) ;
        im = transpose(im) ;
        
    }

}