#include "pch.h"
#include <cmath>
#include "household.h"
#include "parameter.h"
#include <algorithm>

using namespace HOUSEHOLD;
using namespace PARAMETER;
using namespace std;

class household {

const double ALPHA = 0.24;
const double GAMMA = 3.0;

double initial_wealth, working_ability, average_historical_earning, bond_choice ;
int age;

public:
	double maximum_consumption(double labor_choice, double risk_free_return, double equity_return){
		double cash_on_hand;

        return 0;

	}

	double utility(double consumption, double labor) {

		return pow(pow(consumption, ALPHA)*pow(labor, 1.0-ALPHA), 1.0 - GAMMA)/(1.0 - GAMMA);
	}

	double ss_benefit_function() {
		if (age<retire_age) {
			return 0.0;
		}
		else {
			return oasi_adj_factor * (0.9*min(average_historical_earning, zeta1) + 0.32*max(min(average_historical_earning, zeta2) - zeta1, 0.0) + 0.15*max(average_historical_earning - zeta2, 0.0));
		}
	}

};