#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace FEM {
	const static double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899;
	const static double density = 1000;
	const static double YoungModulus = 5e4;
	const static double PoissonRate = 0.49;
	const static double lengthRateLame = YoungModulus / (2 * (1 + PoissonRate));
	const static double volumeRateLame = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
	const static double lengthRate = 4 * lengthRateLame / 3;
	const static double volumeRate = volumeRateLame + 5 * lengthRateLame / 6;
	//static double Hhat = 1e-6;
	//const static double IPC_dt = 0.01;
	//static double Kappa = 0;
}


#endif // !FEM_PARAMETERS_H

