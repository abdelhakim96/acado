

#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>


int main()
{
	USING_NAMESPACE_ACADO


		// INTRODUCE THE VARIABLES:
		// -------------------------
		DifferentialState xB; //Body Position
	DifferentialState xW; //Wheel Position
	DifferentialState vB; //Body Velocity
	DifferentialState vW; //Wheel Velocity

	Control F;
	Disturbance R;

	double mB = 350.0;
	double mW = 50.0;
	double kS = 20000.0;
	double kT = 200000.0;


	// DEFINE A DIFFERENTIAL EQUATION:
	// -------------------------------
	DifferentialEquation f;

	f << dot(xB) == vB;
	f << dot(xW) == vW;
	f << dot(vB) == (-kS * xB + kS * xW + F) / mB;
	f << dot(vW) == (kS * xB - (kT + kS) * xW + kT * R - F) / mW;


	// DEFINE LEAST SQUARE FUNCTION:
	// -----------------------------
	Function h;

	h << xB;
	h << xW;
	h << vB;
	h << vW;
	h << F;

	DMatrix Q(5, 5); // LSQ coefficient matrix
	Q(0, 0) = 10.0;
	Q(1, 1) = 10.0;
	Q(2, 2) = 1.0;
	Q(3, 3) = 1.0;
	Q(4, 4) = 1.0e-8;

	DVector r(5); // Reference
	r.setAll(0.0);


	// DEFINE AN OPTIMAL CONTROL PROBLEM:
	// ----------------------------------
	const double tStart = 0.0;
	const double tEnd = 1.0;

	OCP ocp(tStart, tEnd, 20);

	ocp.minimizeLSQ(Q, h, r);

	ocp.subjectTo(f);

	ocp.subjectTo(-200.0 <= F <= 200.0);
	ocp.subjectTo(R == 0.0);


	// SETTING UP THE REAL-TIME ALGORITHM:
	// -----------------------------------
	RealTimeAlgorithm alg(ocp, 0.025);
	alg.set(MAX_NUM_ITERATIONS, 1);
	alg.set(PLOT_RESOLUTION, MEDIUM);

	GnuplotWindow window;
	window.addSubplot(xB, "Body Position [m]");
	window.addSubplot(xW, "Wheel Position [m]");
	window.addSubplot(vB, "Body Velocity [m/s]");
	window.addSubplot(vW, "Wheel Velocity [m/s]");
	window.addSubplot(F, "Damping Force [N]");
	window.addSubplot(R, "Road Excitation [m]");

	alg << window;


	// SETUP CONTROLLER AND PERFORM A STEP:
	// ------------------------------------
	StaticReferenceTrajectory zeroReference("ref.txt");

	Controller controller(alg, zeroReference);

	DVector y(4);
	y.setZero();
	y(0) = 0.01;

	if (controller.init(0.0, y) != SUCCESSFUL_RETURN)
		exit(1);
	if (controller.step(0.0, y) != SUCCESSFUL_RETURN)
		exit(1);

	return EXIT_SUCCESS;
}