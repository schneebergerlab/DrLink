/* this function analyze/cluster barcodes for regions of a chr to get

   1) distribution of molecule length              - x: molecule_len                      , y: molecule_num
   2) distribution of molecule base coverage       - x: read_num * read_len / molecule_len, y: molecule_num
   3) number of moleculer per partition (optional) - x: molecule_num                      , y: partition_num
   
   input: a 5 column gzipped file converted from bam in a previous step.

*/

bool analyze_molecule(int argc, char* argv[]);
