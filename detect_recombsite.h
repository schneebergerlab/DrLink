/* this functions detects recombination sites in meiosis, based on pooled linked read sequencing 

   inputs: 
      
      1. known markers from both parents
      2. variations called in pool     

   output:

      prediction on recombination sites.
*/

bool detect_recombsite(int argc, char* argv[]);
