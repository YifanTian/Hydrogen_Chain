#include "itensor/all.h"

using namespace itensor;
using namespace std;

int main(int argc, char* argv[]) {

	ofstream output, outfile_N, outfile_correlation;
	output.open("output.txt");
	outfile_correlation.open("correlation.txt");
	outfile_N.open("number_particles.txt");

	//Parse the input file
	if (argc != 2) {
		printfln("Usage: %s inputfile_exthubbard", argv[0]);
		return 0;
	}
	auto input = InputGroup(argv[1], "input");

	auto N = input.getInt("N");
	auto Npart = input.getInt("Npart", N); //number of particles, default is N (half filling)

	auto nsweeps = input.getInt("nsweeps");
	auto t1 = input.getReal("t1", 1);
	auto t2 = input.getReal("t2", 0);
	auto U = input.getReal("U", 0);
	auto V1 = input.getReal("V1", 0);
	auto quiet = input.getYesNo("quiet", false);

	auto table = InputGroup(input, "sweeps");
	auto sweeps = Sweeps(nsweeps, table);
	println(sweeps);

	//
	// Initialize the site degrees of freedom.
	//
	auto sites = Hubbard(N);

	//
	// Create the Hamiltonian using AutoMPO      how to periodic, no interaction
	//
	auto ampo = AutoMPO(sites);
	for (int i = 1; i <= N; ++i) {
		ampo += U, "Nupdn", i;
		//if (i == 10) { ampo += U, "Nupdn", i; }
		//if (i != 11) { ampo += U, "Nupdn", i; }
	}
	
	for (int b = 1; b < N; ++b) {
		ampo += -t1, "Cdagup", b, "Cup", b + 1;
		ampo += -t1, "Cdagup", b + 1, "Cup", b;
		ampo += -t1, "Cdagdn", b, "Cdn", b + 1;
		ampo += -t1, "Cdagdn", b + 1, "Cdn", b;
		//ampo += V1,"Ntot",b,"Ntot",b+1;
	}
	/*
	ampo += -t1, "Cdagup", N, "Cup", 1;
	ampo += -t1, "Cdagup", 1, "Cup", N;
	ampo += -t1, "Cdagdn", N, "Cdn", 1;
	ampo += -t1, "Cdagdn", 1, "Cdn", N;
	*/
	 for(int b = 1; b < N-1; ++b)
	 {
	 ampo += -t2,"Cdagup",b,"Cup",b+2;
	 ampo += -t2,"Cdagup",b+2,"Cup",b;
	 ampo += -t2,"Cdagdn",b,"Cdn",b+2;
	 ampo += -t2,"Cdagdn",b+2,"Cdn",b;
	 }
	 
	auto H = IQMPO(ampo);

	//
	// Set the initial wavefunction matrix product state
	// to be a Neel state.
	//
	auto state = InitState(sites);							//how to 1 particle?
	int p = Npart;
	for (int i = N; i >= 1; --i) {
		if (p > i) {
			println("Doubly occupying site ", i);
			state.set(i, "UpDn");
			p -= 2;
		} else if (p > 0) {
			println("Singly occupying site ", i);
			state.set(i, (i % 2 == 0 ? "Up" : "Dn"));
			//state.set(i, "Up");
			p -= 1;
		} else {
			state.set(i, "Emp");
		}
	}

	auto psi = IQMPS(state);

	Print(totalQN(psi));

	//
	// Begin the DMRG calculation
	//
	auto energy = dmrg(psi, H, sweeps, { "Quiet", quiet });
	
	//
	// Measure correlation function
	//
	int i = 4;
	int j = 8;
	
	for(int i = 1;i<N;i++) {

			j = N;
			//cout<<"i "<<i<<"j "<<j<<endl;
	
			auto adagup_i = sites.op("Adagup",i);
			auto aup_j = sites.op("Aup",j); 
			auto F_i = sites.op("FermiPhase",i);
					
			psi.position(i); 
			
			auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
			auto C = psi.A(i)*F_i;
			//auto mutioperater = F_i*prime(adagup_i);
			C*= prime(adagup_i);
			C*= dag(prime(prime(psi.A(i),Site,2),ir));  //or C = prime(C,Site,-1);
			
			auto correlation = (dag(prime(psi.A(i), Site)) * sites.op("Nup", i) * psi.A(i)).real();
			cout<<"correlation "<<i<<" "<<i<<" "<<correlation<<endl;
			outfile_correlation<<i<<" "<<i<<" "<<correlation<<endl;
		
			auto C1 = C*psi.A(i+1);
			C1*= sites.op("Aup",i+1);
			auto jr = commonIndex(psi.A(i),psi.A(i+1),Link);
			C1*= dag(prime(psi.A(i+1),jr,Site));
			correlation = C1.real();
			cout<<"correlation "<<i<<" "<<i+1<<" "<<correlation<<endl;
			outfile_correlation<<i<<" "<<i+1<<" "<<correlation<<endl;
			//outfile_correlation<<i<<" "<<i+1<<" "<<correlation;
			
			for (int b = i+1; b < j; b++)
			{	
				C*= psi.A(b);
				C*= sites.op("F",b);
				C*= dag(prime(psi.A(b),Link,Site)); 
				
				C1 = C*psi.A(b+1);
				C1*= sites.op("Aup",b+1);
				jr = commonIndex(psi.A(b),psi.A(b+1),Link);
				C1*= dag(prime(psi.A(b+1),jr,Site));
				correlation = C1.real();
				cout<<"correlation "<<i<<" "<<b+1<<" "<<correlation<<endl;
				outfile_correlation<<i<<" "<<b+1<<" "<<correlation<<endl;
				//outfile_correlation<<i<<" "<<b+1<<" "<<correlation;
			}
			//outfile_correlation<<endl;
	}
	psi.position(N);
	auto correlation = (dag(prime(psi.A(N), Site)) * sites.op("Nup", N) * psi.A(N)).real();
	
	cout<<"correlation "<<N<<" "<<N<<" "<<correlation<<endl;
	//outfile_correlation<<N<<" "<<N<<" "<<correlation<<endl;

	Vector upd(N), dnd(N);
	for (int j = 1; j <= N; ++j) {
		psi.position(j);
		upd(j - 1) =
				(dag(prime(psi.A(j), Site)) * sites.op("Nup", j) * psi.A(j)).real();
		dnd(j - 1) =
				(dag(prime(psi.A(j), Site)) * sites.op("Ndn", j) * psi.A(j)).real();
	}
	
	printfln("%d %.10f", N, upd(N-1));
	outfile_correlation<<N<<" "<<N<<" "<<upd(N-1)<<endl;
	outfile_correlation.close();
	
	auto Energy = 0.0;
	// two sites correlation
	for (int b = 1; b < N; b++)
	{		
		psi.position(b); 
		auto ir = commonIndex(psi.A(b),psi.A(b+1),Link);
		auto C = psi.A(b)*sites.op("F",b);
		C*= prime(sites.op("Adagup",b));
		C*= dag(prime(prime(psi.A(b),Site,2),ir));
		
		auto C1 = C*psi.A(b+1);
		C1*= sites.op("Aup",b+1);
		auto jr = commonIndex(psi.A(b),psi.A(b+1),Link);
		C1*= dag(prime(psi.A(b+1),jr,Site));
		auto correlation = C1.real();
		cout<<"two sites correlation"<<b<<" "<<correlation<<endl;
		Energy += correlation;
	}
	cout<<"Energy as sum of correlation"<<" "<<Energy<<endl;
	
	
	Energy = 0.0;
	// three sites correlation
	for (int b = 1; b < N-1; b++)
	{		
		psi.position(b); 
		auto ir = commonIndex(psi.A(b),psi.A(b+1),Link);
		auto C = psi.A(b)*sites.op("F",b);
		C*= prime(sites.op("Adagup",b));
		C*= dag(prime(prime(psi.A(b),Site,2),ir));
		
		C*= psi.A(b+1);
		C*= sites.op("F",b+1);
		C*= dag(prime(psi.A(b+1),Link,Site)); 
		
		auto C1 = C*psi.A(b+2);
		C1*= sites.op("Aup",b+2);
		auto jr = commonIndex(psi.A(b+1),psi.A(b+2),Link);
		C1*= dag(prime(psi.A(b+2),jr,Site));
		auto correlation = C1.real();
		cout<<"two sites correlation"<<b<<" "<<correlation<<endl;
		Energy += correlation;
	}
	cout<<"Energy as sum of correlation"<<" "<<Energy<<endl;
	
	//Entanglement entropy
	psi.position(N/2);
	auto wf = psi.A(N/2)*psi.A((N/2)+1);
	auto S = psi.A(N/2);
	IQTensor V,D;
	auto spectrum = svd(wf,S,V,D);
	
	Real SvN = 0.0;
	
	for(auto p:spectrum.eigs())
	{
		if(p>1E-12) SvN += -p*log(p);
	}
	printfln("Across bond b=%d, SvN = %.10f",N/2,SvN);
	
	//
	// Measure spin densities
	//
	/*
	Vector upd(N), dnd(N);
	for (int j = 1; j <= N; ++j) {
		psi.position(j);
		upd(j - 1) =
				(dag(prime(psi.A(j), Site)) * sites.op("Nup", j) * psi.A(j)).real();
		dnd(j - 1) =
				(dag(prime(psi.A(j), Site)) * sites.op("Ndn", j) * psi.A(j)).real();
	}
	*/

	println("Up Density:");
	for (int j = 0; j < N; ++j)
		printfln("%d %.10f", 1 + j, upd(j));
	println();

	println("Dn Density:");
	for (int j = 0; j < N; ++j)
		printfln("%d %.10f", 1 + j, dnd(j));
	println();

	println("Total Density:");
	for (int j = 0; j < N; ++j)
		printfln("%d %.10f", 1 + j, (upd(j) + dnd(j)));
	println();
	
	for (int j = 0; j < N; ++j)
	{
		outfile_N<<1+j<<" "<<upd(j)<<" "<<dnd(j)<<" "<<(upd(j) + dnd(j))<<endl;
	}
	outfile_N.close();

	//
	// Print the final energy reported by DMRG
	//
	printfln("\nGround State Energy = %.10f, energy per site = %.10f  ", energy, energy/N);

	return 0;
}
