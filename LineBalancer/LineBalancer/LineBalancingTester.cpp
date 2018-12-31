#include "LineBalancingTester.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

double LineBalancingTester::fitness()
{
	double sum = 0.0;
	for (auto& gene : resultPacking)
	{
		int fill = 0;
		for (auto& item : gene)
			fill += items[item];
		sum += std::pow(1.0*fill / binCapacity, 2);
	}
	return sum / resultPacking.size();
}

void LineBalancingTester::displayConsole()
{
	std::cout << "\n";
	for (int i = 0; i < resultPacking.size(); ++i)
	{
		std::cout << "workstation " << i << ": ";
		for (auto& x : resultPacking[i]) std::cout << x << " ";
		std::cout << "\n";
	}
	double result = fitness();
	std::cout << "Result: " << result << "\n";
	std::cout << "Best possible: " << result + 0.005 << "\n";
}

void LineBalancingTester::displayGraphics()
{
	sf::Color blue(63, 123, 178), orange(245, 128, 61), grey(66, 66, 66);

	int binHeight = 5, binGap = 1;
	for (int i = 0; i < resultPacking.size(); ++i)
	{
		sf::RectangleShape workstationRect;
		workstationRect.setSize(sf::Vector2f(std::accumulate(resultPacking[i].begin(), resultPacking[i].end(), 0, [this](int soFar, int curr) {return soFar + items[curr]; }),
											binHeight));
		workstationRect.setFillColor(orange);
		workstationRect.setPosition(0, i * (binHeight + binGap));
		window->draw(workstationRect);
	}
	sf::RectangleShape whiteLine(sf::Vector2f(1, resultPacking.size() * (binHeight+binGap)));
	whiteLine.setFillColor(sf::Color::White);
	whiteLine.setPosition(binCapacity, 0);
	window->draw(whiteLine);

	sf::RectangleShape greyLine(sf::Vector2f(1, resultPacking.size() * (binHeight + binGap)));
	greyLine.setFillColor(blue);
	greyLine.setPosition(capacity, 0);
	window->draw(greyLine);

	if (bestFitness.size() > 1)
	{
		int topY = resultPacking.size() * (binHeight + binGap) + 60;
		int graphHeight = 200, graphWidth = 1000;
		int leftGap = 30;
		sf::RectangleShape oy(sf::Vector2f(1, graphHeight + 5));
		oy.setFillColor(grey);
		oy.setPosition(leftGap, topY);
		window->draw(oy);
		sf::RectangleShape ox(sf::Vector2f(graphWidth + 5, 1));
		ox.setFillColor(grey);
		ox.setPosition(leftGap - 5, topY + graphHeight);
		window->draw(ox);

		double scaleX = graphWidth / (bestFitness.size() - 1), scaleY = graphHeight;
		int oX = leftGap, oY = topY + graphHeight;
		sf::Font font;
		font.loadFromFile("c:\\Users\\mickl\\Documents\\Visual Studio 2017\\Projects\\LineBalancer\\Debug\\consola.ttf");
		for (double d = 0.1; d < 1.008; d += 0.1)
		{
			std::stringstream ss;
			ss << std::fixed << std::setprecision(1) << d;
			sf::Text text(ss.str(), font);
			text.setPosition(2, oY - d*scaleY - 7);
			text.setCharacterSize(12);
			window->draw(text);

			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(oX - 3, oY - d*scaleY), grey),
				sf::Vertex(sf::Vector2f(oX + 3, oY - d*scaleY), grey)
			};
			window->draw(line, 2, sf::Lines);
		}
		for (int i = 0; i < bestFitness.size(); ++i)
		{
			std::stringstream ss;
			ss << i;
			sf::Text text(ss.str(), font);
			text.setPosition(oX + i*scaleX, oY + 2);
			text.setCharacterSize(12);
			window->draw(text);
		}

		for (int i = 0; i < bestFitness.size() - 1; ++i)
		{
			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(oX + i *scaleX, oY - bestFitness[i] * scaleY), orange),
				sf::Vertex(sf::Vector2f(oX + (i + 1)*scaleX, oY - bestFitness[i + 1] * scaleY), orange)
			};
			window->draw(line, 2, sf::Lines);
		}

		//double bestVal = (1.0*capacity*binsAmount) / (binCapacity*resultPacking.size());
		double bestVal = bestFitness.back()+0.005;
		sf::Vertex bestLine[] =
		{
			sf::Vertex(sf::Vector2f(oX,								oY - bestVal*scaleY), blue),
			sf::Vertex(sf::Vector2f(oX + (bestFitness.size() - 1)*scaleX, oY - bestVal*scaleY), blue)
		};
		window->draw(bestLine, 2, sf::Lines);
	}
}

GeneticBalancer::PrecedenceGraph LineBalancingTester::generateAcyclicPrecedenceGraph(int n)
{
	std::vector<std::pair<int, int>> g;
	
	for (int i = 0; i < n-1; ++i)
	{
		int outEdges = (random(0, 99)<70? 1 : 2);
		for (int j = 0; j < outEdges; ++j)
			g.push_back({ i, random(i+1, n-1) });
	}
	if (std::find_if(g.begin(), g.end(), [n](std::pair<int, int> p) { return p.first==n-1 || p.second==n-1; }) == g.end())
		g.push_back({ n-2, n-1 });

	/*
	std::ofstream testFile("c:\\Users\\mickl\\Documents\\Visual Studio 2017\\Projects\\LineBalancer\\Tests\\1.txt");
	//std::cout << "Precedence graph:\n";
	for (auto&x : g) testFile << x.first << " " << x.second << "\n";
	*/
	return GeneticBalancer::PrecedenceGraph(g);
}

void LineBalancingTester::test()
{
	test(random(300, 1000), random(0, 10)); //2-1000, 0-20
}

void LineBalancingTester::test(int binCapacity, double leeway)
{
	capacity = binCapacity * (1.0 - leeway/100.0);
	binsAmount = random(30, 77); //1-100
	std::vector<int> generatedItems;
	for (int i = 0; i < binsAmount; ++i)
	{
		int currentCapacity = capacity;
		while (currentCapacity > 0)
		{
			int curr = random(1, currentCapacity);
			generatedItems.push_back(curr);
			currentCapacity -= curr;
		}
	}
	std::random_shuffle(generatedItems.begin(), generatedItems.end());
	system("CLS");
	std::cout << "Bin capacity: " << binCapacity << "\n"
		<< "Leeway: " << leeway << "\n"
		<< "Capacity: " << capacity << "\n"
		<< "Bins Amount: " << binsAmount << "\n";
	test(generatedItems, binCapacity, generateAcyclicPrecedenceGraph(generatedItems.size()));
}

void LineBalancingTester::test(std::vector<int> items, int binCapacity, const GeneticBalancer::PrecedenceGraph& pg)
{
	//std::cout << "{"; for (auto&x : items) std::cout << x << ","; std::cout << "} ["<<items.size()<<" items]\n";
	this->items = items;
	this->binCapacity = binCapacity;

	/*
	std::ofstream testFile("c:\\Users\\mickl\\Documents\\Visual Studio 2017\\Projects\\LineBalancer\\Tests\\1.txt", std::ios::app);
	testFile << this->binCapacity << "\n";
	for (auto& x : items) testFile << x << " ";
	*/

	long elapsedTime;
	resultPacking = GeneticBalancer().balance(items, binCapacity, pg, bestFitness, elapsedTime);
	
	displayConsole();
	std::cout << "Algorithm running time: " << elapsedTime << "ms\n\n\n\n\n\n\n\n";
	displayGraphics();
}

void LineBalancingTester::run()
{
	srand(time(0));
	test();
}

void LineBalancingTester::run(int i)
{
	std::vector<std::pair<int, int>> readPG;
	std::vector<int> readItems;
	int readBinCapacity;

	std::stringstream fileName;
	fileName << "c:\\Users\\mickl\\Documents\\Visual Studio 2017\\Projects\\LineBalancer\\Tests\\"<<i<<".txt";
	std::ifstream inFile(fileName.str());

	//read from file
	while (true)
	{
		int a, b;
		inFile >> a >> b;
		if (a != -1)
			readPG.push_back({ a, b });
		else
		{
			readBinCapacity = b;
			break;
		}
	}
	for (int a; inFile >> a; readItems.push_back(a));

	system("CLS");
	test(readItems, readBinCapacity, readPG);
}
