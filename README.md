# RRAM_model
A Verilog-A compact RRAM model (designed and developed by Stanford &amp; ASU) which describes the I-V characteristic of this emerging memory device. 

`include "constants.vams"
`include "disciplines.vams"

module RRAM(Nt, Nb);

	inout Nt, Nb;
	electrical Nt, Nb, gap_out, R_out, temp_out;

	real kb =  `P_K;	// Boltzmann's constant= 1.3806503e-23 (J/K)
	real q =  `P_Q;		// Electron charge= 1.6e-19 (C)

	// Device parameters
	parameter real L = 5e-9 from (0:inf);           // Oxide thickness (m)
	parameter real gap_min = 0.1e-9 from (0:L);		// Min. gap distance (m)
	parameter real gap_max = 1.7e-9 from (gap_min:L);	// Max. gap distance (m)
	parameter real gap_ini = 0.1e-9 from [gap_min:gap_max];	// Initial gap distance (m)
	parameter real a0 = 0.25e-9 from (0:inf);		// Atomic distance (m)
	parameter real Eag = 1.501 from (0:inf);  	// Activation energy for vacancy generation (eV)
	parameter real Ear = 1.5 from (0:inf);  	// Activation energy for vacancy recombination (eV)

	// I-V characteristics
	parameter real I0 = 6.14e-5 from (0:inf);
	parameter real g0 = 2.7505e-10 from (0:inf);	
	parameter real V0 = 0.43 from (0:inf);

	// Gap dynamics
	parameter real Vel0 = 150 from (0:inf); 
	parameter real gamma0 = 16.5 from (0:inf);
	parameter real g1 = 1e-9 from (0:inf);
	parameter real beta = 1.25 from (0:gamma0/(pow(gap_max/g1,3)));
	real gap;       // Current gap (m)
	real gap_ddt;   // Gap growth velocity (m/s)
	real gamma;     // Local enhancement factor

	// Temperature dynamics
	parameter real T0 = 273+25 from (0:inf);		// Ambient temperature (K)
	parameter real Cth = 3.1825e-16 from (0:inf);	// Effective thermal capacitance (J/K)
	parameter real Tau_th = 2.3e-10 from (0:inf);	// Effective thermal time constant (s)
	real temperature;   // Current temperature (K)

	// Simulation time control
	parameter real tstep = 1e-9 from (0:inf);	// Max. internal timestep (s)
	real c_time;	// Current timestep (s)
	real p_time;	// Previous timestep (s)
	real dt;		// Difference between current and previous timestep (s)
	
	// Resitance
	parameter real Vread = 0.1; // Read voltage (V)
	real Iread;     // Read current @Vread (A)


	analog begin

		@(initial_step)  begin
			temperature = T0;
			gap = gap_ini;
		end

		$bound_step(tstep);	// Bound the time step
		c_time = $abstime;	// Get current timestep
		dt = c_time-p_time;
		
		// Calculate the local enhancement factor
		gamma = gamma0 - beta*pow(gap/g1, 3);	

		// Gap dynamics
		gap_ddt = -Vel0*( exp(-q*Eag/kb/temperature)*exp(gamma*a0/L*q*V(Nt,Nb)/kb/temperature) - exp(-q*Ear/kb/temperature)*exp(-gamma*a0/L*q*V(Nt,Nb)/kb/temperature) );		
		gap = gap + gap_ddt*dt;
		
		// Limitation on the gap growth
		if (gap > gap_max)
			gap = gap_max;
		else if (gap < gap_min)
			gap = gap_min;
		
		// Calculate the current
		I(Nt,Nb) <+ I0*exp(-gap/g0)*sinh(V(Nt,Nb)/V0);

		// Calculate the local temperature (implicit form)
		temperature = (temperature + dt*(abs(V(Nt,Nb)*I(Nt,Nb))/Cth+T0/Tau_th))/(1+dt/Tau_th);

		// Parameter extraction
		V(gap_out) <+ gap;			// Gap distance
		V(temp_out) <+ temperature;	// Temperature
		Iread = I0*exp(-gap/g0)*sinh(Vread/V0);
		V(R_out) <+ Vread/Iread;	// Resistance @Vread
		
		p_time = $abstime;	// Record current timestep for the previous timestep at next timestep

	end
endmodule
