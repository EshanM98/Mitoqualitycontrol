This folder contains the code for the quantitative models of the project and an excel sheet of all the experimental data 

functions.py contains all the functions used to in analyzing the degradation and replication models

	The binning function is used to split the experimental data into some number of bins and then average each bin with a minimum number of data points
		Inputs:
			bins - the number of bins to split experimental data
			binthreshold - the minimum number of data points required for a bin to be averaged
			prevgen_all - the data points for the mutant genome in the previous generation to be binned and averaged
			nextgen_all - the data points for the mutant genome in the next generation to be binned and averaged
		Output:
			groupsprev_avg - The average mutant genome in the previous generation for each bin that has a number of datapoints greater than or equal to the binthreshold
			groupsnext_avg - The average mutant genome in the next generation for each bin that has a number of datapoints greater than or equal to the binthreshold
			groups_prev - The data points for the mutant genome in the previous generation found in each bin
			groups_next - The data points for the mutant genome in the next generation found in each bin


	deg2A is used to predict the mutant genome in the next generation via the unmodified degradation model in region 2A
		Input: 
			f0 - float value of the mutant genome in the previous generation
		Output:
			f_next - float value of the mutant genome in the next generation
			t2 - Total number of mtDNA copies 
			t1 - Total number of mutant mtDNA copies isolated in mitochondria with only mutant mtDNA

	deg1 is used to predict the mutant genome in the next generation via the unmodified degradation model in region 1
		Input: 
			f0 - float value of the mutant genome in the previous generation
		Output:
			f_next - float value of the mutant genome in the next generation
			t2 - Total number of mtDNA copies 
			t1 - Total number of mutant mtDNA copies isolated in mitochondria with only mutant mtDNA

	inter_reg is used to build an interpolated model between the unmodified degradation model of region 1 and 2A
		Input:
			f0 - float value of mutant genome in the previous generation
			p - interpolation value
		Output:
			f_inter - interpolated float value of mutant genome in the next generation

	SSE is used to determine the sum of squared errors for the unmodified degradation model and interpolated model
		Input:
			prevgen - Mutant genome in the previous generation (using binned and averaged experimental data)
			nextgen - Mutant genome in the next generation (using binned and averaged experimental data)
			p - interpolation value used for interpolated degradation model
		Output:	
			sse_r1 - sum of squared errors value for unmodified region 1 degradation model
			sse_r2 - SSE value for unmodified region 2a degradation model
			sse_int - SSE value for interpolated degradation model

	deg2A_onerem and deg2A_tworem predict the modified degradation model for the removal of mutant mitochondria with only one nucleoid, or mutant mitochondria with only one or two nucleoids
		Inputs: 
			f0 - mutant genome in the previous generation
		Output:
			f_next - predicted mutant genome in the next generation

	deg2A_int determines the interpolation between the modified region 2A degradation model between selection of mutant mitochondria with one nucleoid and mutant mitochondria with one or two nucleoids
		Inputs:
			f0 - mutant genome in the previous generation
			p - interpolation value
		Output:
			f_inter - predicted mutant genome in the next generation of interpolated model

	deg1_onerem, deg1_tworem and deg1_threerem predict the modified degradation model for the removal of mutant mitochondria with only one nucleoid, mutant mitochondria with only one or two nucleoids, or mutant mitochondria with one, two or three nucleoids
		Inputs: 
			f0 - mutant genome in the previous generation
		Output:
			f_next - predicted mutant genome in the next generation

	rep_comp returns the solution to the complementation replication model (Eq 2.7 of thesis)
		Input:
			f1 - mutant genome isolated to mitochondria with only mutant mtDNA
			f2 - all mutant genome 
			n - number of replication cycles
			lamb - replication factor for mutant mtDNA
	rep_noncomp returns the solution to the no complementation replication model (Eq. 2.9 of thesis)
		Input:
			f - mutant genome in the previous generation
			n - number of replication cycles
			lamb - replication factor for mutant mtDNA

	rep_nocomp is used to predict the mutant genome in the next generation following the no complementation replication model
		Input:
			f0 - mutant genome in the previous generation
			lamb - replication factor for mutant mtDNA
			mtc_fin - final count of mtDNA reached from replication 
		Output:
			f_nextgen_nocomp - mutant genome in the next generation
	rep_r1 is used to predict the mutant genome in the next generation following the region 1 replication model with complementation
		Input:
			f0 - mutant genome in the previous generation
			lamb - replication factor for mutant mtDNA
			mtc_fin - final count of mtDNA reached from replication 
		Output:
			f_nextgen_comp- mutant genome in the next generation
			n - number of replication cycles needed to reach final mtDNA count

	rep_r2 is used to predict the mutant genome in the next generation following the region 2a replication model with complementation
		Input:
			f0 - mutant genome in the previous generation
			lamb - replication factor for mutant mtDNA
			mtc_fin - final count of mtDNA reached from replication 
		Output:
			f_nextgen_comp- mutant genome in the next generation
			n - number of replication cycles needed to reach final mtDNA count

	sse_rep_nocomp, sse_rep_r1 and sse_rep_r2 determine the SSE for each replication model
		Inputs:
			lamb - replication factor for mutant mtDNA
			prevgen - mutant genome in the previous generation (from binned data)
			nextgen - mutant genome in the next generation (from binned data)
			mtc_fin - final count of mtDNA needed from replication model
		Outputs:
			sse_noncomp - SSE of no complementation replication model
			sse_comp - SSE of complementation replication models

Experimental data.py is used to plot all experimental data as shown in Fig 2.1 of thesis

Degradation Model - Unmodified.py is used to plot the unmodified degradation model (Fig 3.1A of thesis)
	bins and binthresh can be changed to alter the binning of experimental data
	p can be changed to alter the interpolation between region 1 and 2a

Degradation Model - limited removal.py is used to plot the selective degradation model (Fig 3.1B of thesis)
	bins and binthresh can be changed to alter the binning of experimental data
	p can be changed to alter the interpolation between region 1 and 2a
	p_int can be changed to alter the interpolation between the region 2a selective degradation model
	
	The unmodified degradation model is loaded from an excel file (as deg_model) and can be plotted to compare to the selective model

Interpolated degradation model.py is used to determine the minimum SSE for the interpolation model of the unmodified degradation model and the selective degradation model

Replication model.py is used to plot the replication model (Fig 3.2C of thesis)
	bins and binthresh can be changed to alter the binning of experimental data
	lamb_nocomp, lamb_r1 and lamb_r2 indicate the best replication factors for mutant mtDNA for up to six replication cycles
	mtc_fin indicates the final number of mtDNA copies for each replication factor

Minimizing lambda for replication model.py is used to plot the minimized SSE and lambda values for each number of replications (Fig 3.2 A and B of thesis)

Final models.py is used to plot the best quantitative models against experimental data (Fig 3.3 of thesis)
