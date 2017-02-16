//Example of how to use the Logger

#include "Logger.hh"

int main() {

	try {

		//Tell the logger to only show message at INFO level or above
		Logger::Instance()->setLogLevel(Logger::INFO); 
	
		//Tell the logger to throw exceptions when ERROR messages are received
		Logger::Instance()->enableCriticalErrorThrow();

		//Send an INFO message (messages are just strings)
		Logger::Instance()->write(Logger::INFO,"Hello from Logger");

		//Send an INFO message with a number in it (make the string using a stringstream)
		std::stringstream msg;
		msg << "5 + 5 = " << (5+5);
		Logger::Instance()->write(Logger::INFO,msg.str());

		//Another way to send an INFO message with a number, using std::to_string to turn a number into a string
		Logger::Instance()->write(Logger::INFO,"1 + 6 = " + std::to_string(1+6));

		//Send an ERROR message, this will terminate the program
		Logger::Instance()->write(Logger::ERROR,"Something terrible happened");

	}

	//Catch Logger exceptions
	catch (CriticalError& e) {
		std::cerr << "A critical error occurred, exiting" << std::endl;
		return -1; //Exit program wth an error code
	}
 
	return 0;

}