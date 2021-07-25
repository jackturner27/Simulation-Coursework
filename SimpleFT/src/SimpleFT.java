public class SimpleFT {

	public static int N = 256;

	public static void main(String [] args) throws Exception {
  
          double [] [] X = new double [N] [N] ;
          ReadPGM.read(X, "wolf.pgm", N) ;
  
          DisplayDensity display =
                  new DisplayDensity(X, N, "Original Image") ;
          
          double starttime = System.currentTimeMillis();
  
          double [] [] CRe = new double [N] [N], CIm = new double [N] [N] ;
  
          for(int k = 0 ; k < N ; k++) {
              for(int l = 0 ; l < N ; l++) {
                  double sumRe = 0, sumIm = 0 ;
                  // Nested for loops performing sum over X elements
                  for(int m = 0; m < N; m++) {
                      for(int n = 0; n < N; n++) {
                           double arg = (-2 * Math.PI * ((m * k) + (n * l))) / N ;
                           double cos = Math.cos(arg) ;
                           double sin = Math.sin(arg) ;
                           sumRe += cos * X [m] [n] ;
                           sumIm += sin * X [m] [n] ;
                      }
                  }
                  CRe [k] [l] = sumRe ;
                  CIm [k] [l] = sumIm ;
              }
              System.out.println("Completed FT line " + k + " out of " + N) ;
          }
  
          Display2dFT display2 =
                  new Display2dFT(CRe, CIm, N, "Discrete FT") ;
      
          
          
          int cutoff = N/8 ;  // for example
          for(int k = 0 ; k < N ; k++) {
              int kSigned = k <= N/2 ? k : k - N ;
              for(int l = 0 ; l < N ; l++) {
                  int lSigned = l <= N/2 ? l : l - N ;
                  if(Math.abs(kSigned) < cutoff || Math.abs(lSigned) < cutoff) {
                      CRe [k] [l] = 0 ;
                      CIm [k] [l] = 0 ;
                  }
              }
          }

          Display2dFT display2a =
                  new Display2dFT(CRe, CIm, N, "Truncated FT") ;
          
         
          double [] [] reconstructed = new double [N] [N] ;

          for(int m = 0 ; m < N ; m++) {
              for(int n = 0 ; n < N ; n++) {
                  double sum = 0;
                  //... nested for loops performing sum over C elements
                  for(int k = 0; k < N; k++) {
                	  for(int l = 0; l < N; l++) {
                		  
                		  double arg = (2 * Math.PI * ((k * m) + (n * l))) / N ;
                		  double cos = Math.cos(arg);
                		  double sin = Math.sin(arg);
                		  
                		  sum += (cos * CRe[k][l]) - (sin * CIm[k][l]);
                	  }
                  }
                  reconstructed [m] [n] = sum ;
              }
              System.out.println("Completed inverse FT line " + m + " out of " + N) ;
          }

          DisplayDensity display3 =
                  new DisplayDensity(reconstructed, N, "Reconstructed Image") ;
          
          double endtime = System.currentTimeMillis();
          System.out.println(endtime - starttime + "ms");
	}
	
}