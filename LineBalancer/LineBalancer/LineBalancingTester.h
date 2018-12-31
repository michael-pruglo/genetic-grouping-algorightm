#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include "GeneticBalancer.h"

class LineBalancingTester
{
	std::vector<int> items;
	int binCapacity, capacity, binsAmount;
	std::vector<std::vector<int>> resultPacking;
	std::vector<double> bestFitness;

	int random(int min, int max) { return rand() % (max - min + 1) + min; }
	double fitness();
	void displayConsole();
	void displayGraphics();

	GeneticBalancer::PrecedenceGraph generateAcyclicPrecedenceGraph(int n);
	void test();
	void test(int binCapacity, double leeway);
	void test(std::vector<int> items, int binCapacity, const GeneticBalancer::PrecedenceGraph& pg);
public:
	void run();
	void run(int i);

private:
	sf::RenderWindow* window;
public:
	LineBalancingTester(sf::RenderWindow* w) : window(w) {};
};

