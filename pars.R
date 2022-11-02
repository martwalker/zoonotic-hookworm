# parameters
pars = c(R0=4, #basic reproduction number
         omega11=0.75, #relative intensity of intra-host transmission R011/R022 
         omega21=0.75, #relative intensity of inter-host transmission R021/R012
         muW=1/1.5,  #mortality rate adult hookworms
         muL=1/0.5,  #mortality rate of hookworm larvae (depricated)
         b=0.25,     #strength of density dependent transmission
         mu1=0.2,    #mortality rate of host1
         mu2=0.02,   #mortality rate host 2
         a = 420.7,  #maximum fecundity of female hookworms in absence of density-dependent constraints
         z=2000,     #epg threshold of moderate intensity infection (not used)
         rho=0.5,    #proportion of hookworms that are female
         k1=.25,     #overdispersion of hookworms among host 1
         k2=0.25,    #overdispersion of hookworms among host 2
         start.tx=100,#treatment start time
         epsilon1 = 1,#efficacy of treatment in host 1
         epsilon2=0.9,#efficacy of treatment in host 2
         c1 = 0.5,    #treatment coverage in host 1
         c2=0.3,     #treatment coverage in host 2
         freq.tx1=1, #frequency of treatment in host 1 (n per year = 1/freq)
         freq.tx2=1, #frequency of treatment in host 2 (n per year = 1/freq)
         n.tx1=2,    #number of treatments in host 1
         n.tx2=2,    #number of treatments in host 2
         stop.t=120, #stop time
         dt=1/12,    #step size
         dtout = 12/12, #output step size
         equib=1,    #toggle equilibrium or interventions
         tol=1.0E-6, #internal control parameter
         block=100,
         niter=100)
