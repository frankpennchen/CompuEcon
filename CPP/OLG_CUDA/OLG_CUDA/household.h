#pragma once

namespace HOUSEHOLD
{

	class household {

	double initial_wealth, working_ability, average_historical_earning, bond_choice;
	int age;

	public:
		double maximum_consumption(double labor_choice, double risk_free_return, double equity_return);
		double utility(double consumption, double labor);
		double ss_benefit_function();
	};
}

