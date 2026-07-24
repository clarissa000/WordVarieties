# Summer Research 2026

Project on word varieties of one-relator groups.

To run any code that imports sage, you must use sage file_name, rather than python3

The main files to look at are database_dim_func.py which computes the csv (including the Groebner basis) and dimension_random_func.py which computes the dimension of the Zariski-dense variety for a sample of random words of a specified length

The database itself for words up to length 16 is databases/word_dimensions_nogb.csv. This stored the computed dimension of the variety defined by the ideal <pw -2, p(aw) - x, p(bw) - y>, the dimension of the Zariski-dense part and the cardinality if the variety is finite (both with and without multiplicity). It is not all word - but each unique variety using the function canon.py to convert to a canonical form. 

All files in analysis are for analysing this database and looking at dimension 2 words.

All files in linton are the adjustments to add this data to the database created by Marco Linton, accessable here: https://warwick.ac.uk/fac/sci/maths/people/staff/linton/homepage/ for all with two generators

All files in trace_and_dimension_calculations are different algorithms to compute these. 