#pragma once

namespace HOUSEHOLD
{

	class household {

	const double ALPHA = 0.24;
	const double GAMMA = 3.0;

	double initial_wealth, working_ability, average_historical_earning, bond_choice;
	int age;

	public:
		double maximum_consumption(double labor_choice, double risk_free_return, double equity_return);
		double utility(double consumption, double labor);

	};
}

