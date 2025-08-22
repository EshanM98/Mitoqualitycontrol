This folder contains all of the python scripts for simulating piecemeal mitophagy

func_gillespie_stochastic.py contains the functions used to run the simulation
	mito_func is used to determine the active mitophagy rate by solving Eq. 2.23a and 2.23b of thesis
		Inputs: 
			kp - production rate
			kd_p - prime rate given by Eq. 2.23b
			f - bite size
			t - time (set as 2)
			N0 - initial protein count (set as 10000)
		Outputs:
			Nt - protein count after time t
			km_active - active mitophagy rate

	steady is used to determine the decay rate by solving Eq. 2.22 from the thesis
		Inputs:
			km_back - background mitophagy rate
			f - bite size
		Outputs:
			kd - decay rate

	gillespie_stoch used to track protein counts in a mitochondrial fragment
		Inputs:
			protein_init - initial protein count in a mito fragment (set as 10000 proteins)
			mut_init - Initial fraction of mutant proteins in a mitochondrial fragment (set as 0.5)
			prodrate - production rate (set as protein_init/19)
			mtc_count - number of mitochondrial fragments (set as 1000)
			end_time - time at which Gillespie algorithm ends
			tau - time of simulation (same time as t in mito_func)
			protein_split - mtDNA presence: 0 means only wtDNA prod, 1 means only mutant prod
			mito - Active mitophagy variable: 0 means active mitophagy is off, 1 means its on
			back_mito - background mitophagy rate (set as 0.01)
			threshold - Minimum wildtype protein fraction at which simulation ends (Fraction of current wildtype/Initial protein count)
			active_prod - fraction of prodrate during active mitophagy
			drop - factor by which protein_init drops in tau
			size_dist - array of mitochondrial sizes
			prob_dist - array of probabilities corresponding to mitochondrial sizes
			bite_stdev - standard deviation of bite size distribution
		Outputs:
			mtcm_wcount - dataframe of wildtype protein counts
			mtcm_mcount - dataframe of mutant protein counts
			mtcm_totcount - dataframe of total protein counts
			tm_tot - dataframe of timesteps of simulation
			t_fin - end time of simulation
			mito_rec - list of recorded mitophagy rates for each mitochondria
			decay_rec - list of recorded decay rates for each mitochondria
			bite_size_rec - average bite size of every mitophagy event in a mitochondrial fragment
			bite_draw_rec - bite size for every mitochondrial fragment

prob_func.py contains functions called in various codes to plot data
	
	size_lim_eq - Distribute probability under min size equally throughout all sizes about min size 
		Inputs:
			size_prob - array of the probabilities of mitochondrial sizes (gaussian distribution)
			f_size - all sizes of mitochondria
			f_lim - smallest size of mitochondria
		Outputs:
			p_rec - array of probabilities 
			updated_size - array of sizes of mitochondria

	init - initializes and normalizes size distribution (uses size_lim_eq)
		Inputs:
			x - all sizes of mitochondria
			prob - probability of size distribution
			mu - mean of gaussian function
			sigma - standard deviation of gaussian function
			f_lim - minimum size of mitochondrial fragment
		Outputs:
			gauss_fix - array of probability of initialized distribution
			size_fix - array of sizes of initialized distribution
			y_int - cumulative integration of distribution
	
	choose_size - selects a size from initialized size distribution
		Inputs
			size - array of sizes from initialized distribution
			prob - array of probabilities of initialized distribution
		Outputs:
			size-choose - chosen size

	Poisson - return solution to poisson function
		Input: 
			lamb - expected number of events
			k - number of piecemeal mitophagy events

	Norm - return solution to normal function
		Input:
			lamb - expected number of events
			k - number of piecemeal mitophagy events

Prob fun deterministic.py is used to plot the deterministic probability distributions (Fig 3.4A)

	bite can be varied to examine changes in the probability distribution

Size tests.py is used to plot the size distribution (Fig 3.4B)
	
	The size of the distribution (size), minimum mitochondria size (f_lim), range of sizes (x) and standard deviation (sigma) can be altered)

Single probability distribution.py is used to plot the stochastic probability distribution of mitochondrial fragments (Fig 3.4C)
	
	The variables to initialize the size distributions can be altered (variables mentioned in Size tests.py) and the inputs to the gillespie_stoch can be altered

Threshold plot.py is used to plot Fig 3.5A, determining the probability distributions for various bite sizes
	The same variables can be altered as in Single probability distribution.py

prob dip for threshold plots.py is used to plot Fig 3.5B and C

heat map.py is used to plot Fig 3.6AB
	Fig 3.6CD were plotted by iterating through drops instead of active_prod
	Set drop = 0.5, active_prod = np.linspace(0,1,10), change line 66 to "for m in range(len(active_prod))", swap indices from drop[m] to active_prod[m] 
		e.g. line 74
			mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin1,mitorate,decay,bite_size_true,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop[m],size_dist[0],prob_dist[0],std_bite)

			Changed to

			mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin1,mitorate,decay,bite_size_true,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod[m],drop,size_dist[0],prob_dist[0],std_bite)
			 
	The same variables can be altered as in Single probability distribution.py

dimensionless.py is used to plot Fig 3.7


			
	
