#pragma once
#ifndef _EIGEN_DATA_H_
#define _EIGEN_DATA_H_

namespace __GEIGEN__ {

	struct Vector9 {
		double v[9];
	};

	struct Vector12 {
		double v[12];
	};

	struct Matrix3x3d {
		double m[3][3];
	};

	struct Matrix9x9d {
		double m[9][9];
	};

	struct Matrix12x12d {
		double m[12][12];
	};

	struct Matrix9x12d {
		double m[9][12];
	};

	struct Matrix12x9d {
		double m[12][9];
	};
}


#endif // !_EIGEN_DATA_H_

