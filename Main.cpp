#include "pch.h"
#include <iostream>
#include <cmath>
#include <locale>
#include <fstream>
#include "Vectors.h"

using namespace std;

int main(){
	setlocale(LC_ALL, "bulgarian");
	
	int operationtype;
	cout << "Изберете начин за работа с програмата: 1 - от конзолата, 2 - от файл." << endl;
	cin >> operationtype;


	if (operationtype == 1) {
		ConsoleMenu();
	}

	else if (operationtype == 2) {
		FileRead();
	}

	return 0;
}
