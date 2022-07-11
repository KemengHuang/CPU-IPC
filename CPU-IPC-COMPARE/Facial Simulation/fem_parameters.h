#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace FEM {
	const static double PI = 3.1415926535897932;
	const static double density = 1000;
	const static double YoungModulus = 1e4;
	const static double PoissonRate = 0.45;
	const static double explicit_time_step = 0.000001;
	const static double implicit_time_step = 0.0001;
	const static double lengthRateLame = YoungModulus / (2 * (1 + PoissonRate));
	const static double volumeRateLame = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
	const static double lengthRate = 4 * lengthRateLame / 3;
	const static double volumeRate = volumeRateLame + 5 * lengthRateLame / 6;
	const static double friction = 0.2;
	//static double Hhat = 1e-6;
	const static double IPC_dt = 0.01;
	//static double Kappa = 0;
}


#endif // !FEM_PARAMETERS_H

