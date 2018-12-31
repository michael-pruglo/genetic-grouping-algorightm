#include "GeneticBalancer.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <ctime>

GeneticBalancer::PrecedenceGraph::PrecedenceGraph(const std::vector<std::pair<int, int>>& edges)
{
	std::cout << "Building precedence matrix...";
	int n = -1;
	for (auto& e : edges)
	{
		if (n < e.first) n = e.first;
		if (n < e.second) n = e.second;
	}
	++n;

	longestPath.resize(n);
	for (auto& v : longestPath) v.resize(n);

	for (auto& e : edges)
	{
		int u = e.first, v = e.second;
		if (u > v) std::swap(u, v);
		longestPath[u][v] = 1;
	}

	for (int i = n-2; i >= 0; --i)
		for (int j = n - 1; j > i; --j)
			for (int k = i + 1; k < j; ++k)
				if (longestPath[i][k] && longestPath[k][j])
					longestPath[i][j] = std::max(longestPath[i][j], longestPath[i][k] + longestPath[k][j]);

	std::cout << "done.\n";
	/*std::cout << "\nPrecedence graph longest path matrix:\n";
	std::cout << "    "; for (int i = 0; i < n; ++i) std::cout << std::setw(3) << i; std::cout << "\n";
	std::cout << "    "; for (int i = 0; i < n; ++i) std::cout << "---"; std::cout << "\n";
	for (int i = 0; i < n; ++i)
	{
		std::cout << std::setw(3) << i << "|";
		for (int j = 0; j < n; ++j)
			std::cout << std::setw(3) <<longestPath[i][j];
		std::cout << "\n";
	}*/
}

GeneticBalancer::Chromosome GeneticBalancer::randomChromosome()
{
	Chromosome result(*this);
	result.genes.push_back({});

	std::vector<int> itemsIndexes(items.size());
	for (int i = 0; i < itemsIndexes.size(); ++i)
		itemsIndexes[i] = i;
	std::random_shuffle(itemsIndexes.begin(), itemsIndexes.end());
	int currentBin = 0, currentFill = 0;
	for (int i = 0; i < itemsIndexes.size(); )
	{
		int currentItemSize = items[itemsIndexes[i]];
		if (currentFill + currentItemSize <= binCapacity)
		{
			result.genes[currentBin].push_back(itemsIndexes[i]);
			currentFill += currentItemSize;
			++i;
		}
		else
		{
			++currentBin;
			result.genes.push_back({});
			currentFill = 0;
		}
	}
	result.calcFitness();

	return result;
}

double GeneticBalancer::randomZeroToOne()
{
	return rand() / (RAND_MAX + 1.);
}

int GeneticBalancer::random(int min, int max)
{
	return rand() % (max - min + 1) + min;
}

int GeneticBalancer::spinRoulette(const std::vector<double>& probabilities)
{
	int i = 0;
	for (double r = randomZeroToOne(); r > probabilities[i]; ++i);
	return i;
}

/**
Generates initial population
The population is generated randomly and then sorted according to fitness
The first individual is the fittest
*/
std::vector<GeneticBalancer::Chromosome> GeneticBalancer::initPopulation(int size)
{
	std::vector<Chromosome> result;
	result.reserve(size);
	for (int i = 0; i < size; ++i)
		result.push_back(randomChromosome());
	sortPopulation(result);
	return result;
}

/**
Roulette Wheel Selection
Population is guaranteed to be sorted
The first individual is the fittest
*/
std::pair<GeneticBalancer::Chromosome, GeneticBalancer::Chromosome> GeneticBalancer::selectParents(const std::vector<Chromosome>& population)
{
	//prepare the Roulette
	std::vector<double> probabilities(population.size());
	double sum = std::accumulate(population.begin(), population.end(), 0.0, [](double sumSoFar, Chromosome current) {return sumSoFar + current.getFitness(); });
	for (int i = 0; i < population.size(); ++i)
		probabilities[i] = population[i].getFitness()/sum + (i?probabilities[i-1]:0);
	probabilities.back() = 1.0;

	//spin the Roulette
	int id1 = spinRoulette(probabilities), id2;
	do
	{
		id2 = spinRoulette(probabilities);
	} while(id1 == id2);

	return {population[id1], population[id2]};
}

/**
Stochastic universal sampling
*/
std::vector<GeneticBalancer::Chromosome> GeneticBalancer::selectParents2(const std::vector<Chromosome>& population)
{
	double sum = std::accumulate(population.begin(), population.end(), 0.0, [](double sumSoFar, Chromosome current) {return sumSoFar + current.getFitness(); });
	double dist = sum / (POPULATION_SIZE*CROSSOVER_RATE);
	double start = randomZeroToOne() * dist;

	std::vector<Chromosome> keep;
	std::vector<double> fitnessSum(population.size());
	for (int i = 0; i < population.size(); ++i)
		fitnessSum[i] = population[i].getFitness() + (i?fitnessSum[i-1]:0);
	std::cout << "Kept ";
	for (int i = 0; i < (POPULATION_SIZE*CROSSOVER_RATE); ++i)
	{
		double p = start + i*dist;
		int j;
		for (j = 0; fitnessSum[j] < p; ++j);
		keep.push_back(population[j]);
		std::cout << j << " ";
	}
	std::cout << "\n";
	std::random_shuffle(keep.begin(), keep.end());
	return keep;
}

void GeneticBalancer::firstFit(const std::vector<int>& candidates, Chromosome & chromosome)
{
	for (auto& item : candidates)
	{
		bool emplacedFlag = false;
		for (int i = 0; i < chromosome.genes.size(); ++i)
			if (std::accumulate(chromosome.genes[i].begin(), chromosome.genes[i].end(), 0, [this](int sumSoFar, int i) {return sumSoFar + items[i]; }) + items[item] < binCapacity)
			{
				bool cycleDangerFlag = false;
				for (auto& alreadyPacked : chromosome.genes[i])
					if (precedenceGraph.getMaxDistance(item, alreadyPacked) > 1)
						cycleDangerFlag = true;
				if (cycleDangerFlag)
				{
					//std::cout << "--------------------------------------Cycle danger: " << item << " in [";
					//for (auto& alreadyPacked : chromosome.genes[i]) std::cout << alreadyPacked << ",";
					//std::cout << "]\n";

					continue;
				}

				chromosome.genes[i].push_back(item);
				emplacedFlag = true;
				break;
			}
		if (!emplacedFlag)
			chromosome.genes.push_back({ item });
	}
}

void GeneticBalancer::sortPopulation(std::vector<Chromosome>& population)
{
	std::sort(population.begin(), population.end(), [](Chromosome a, Chromosome b) { return a.isFitter(b); });
	
	//std::sort(population.begin(), population.end(), [](Chromosome a, Chromosome b) { return b.isFitter(a); });
}

void GeneticBalancer::printPopulation(const std::vector<Chromosome>& population, int id)
{
	std::cout << "\n\n\n\n\n\n------ " << id << " ---------------------------------------------------------------------------------------------------------------------------------------------\n";
	int maxLen = 0;
	for (auto& item : population)
	{
		if (maxLen < item.toString().size())
			maxLen = item.toString().size();
	}
	for (auto& item : population)
		std::cout << item.toString() << std::string(maxLen-item.toString().size()+2, ' ') << item.getFitness() << "\n";
	std::cout << "\n";
}

void GeneticBalancer::Chromosome::calcFitness()
{
	double sum = 0.0;
	for (auto& gene : genes)
	{
		int fill = 0;
		for (auto& item : gene)
			fill += parent.items[item];
		sum += std::pow(1.0*fill / parent.binCapacity, 2);
	}
	fitness = sum / genes.size();
}

GeneticBalancer::Chromosome & GeneticBalancer::Chromosome::operator=(const Chromosome & b)
{
	//TODO: if (parent != b.parent) throw exception;
	genes = b.genes;
	fitness = b.fitness;
	return *this;
}

bool GeneticBalancer::Chromosome::isFitter(const Chromosome & b) const
{
	return this->fitness > b.fitness;
}

bool GeneticBalancer::Chromosome::isMaximallyFit() const
{
	return fitness == 1.0;
}

GeneticBalancer::Chromosome GeneticBalancer::crossover(const GeneticBalancer::Chromosome & parent1, const GeneticBalancer::Chromosome & parent2)
{
	Chromosome child = parent2;
	
	///1.Select Crossing Section
	int left = random(0, parent1.genes.size()-1), right = random(0, parent1.genes.size() - 1);
	if (left > right) std::swap(left, right);

	///3.Eliminate Doubles & 4.1.Identify Affected Items
	std::vector<int> affectedItems;
	for (int i = left; i <= right; ++i)
	{
		for (auto& item : parent1.genes[i])
		{
			bool foundFlag = false;
			for (int j = 0; foundFlag==false && j<child.genes.size(); ++j)
				for (int k = 0; k < child.genes[j].size(); ++k)
					if (child.genes[j][k] == item)
					{
						affectedItems.insert(affectedItems.end(), child.genes[j].begin(), child.genes[j].end());
						child.genes.erase(child.genes.begin() + j);
						foundFlag = true;
						break;
					}

			affectedItems.erase(std::remove_if(
				affectedItems.begin(), affectedItems.end(),
				[item](int x) { return x==item; }
			), affectedItems.end());
		}
	}

	///2.Insert Crossing Section Groups
	child.genes.insert(child.genes.begin() + random(0, child.genes.size()),
		parent1.genes.begin() + left,
		parent1.genes.begin() + right + 1);
	
	///4.2.Redistribute Affected Items [using FFD]
	std::sort(affectedItems.begin(), affectedItems.end(), [this](int a, int b) { return items[a]>items[b]; });
	firstFit(affectedItems, child);

	///recalculate fitness
	child.calcFitness();
	/*
	double pf = std::max(parent1.getFitness(), parent2.getFitness());
	if (child.getFitness() - pf > 0.1 * (1.0 - pf))
	{
		std::cout << "CROSSOVER: left=" << left << " right=" << right << " affected items ";
		for (auto& x : affectedItems) std::cout << x << " ";
		std::cout<<"\n"
			<< parent1.toString() << " " << parent1.getFitness() << "\n" 
			<< parent2.toString() << " " << parent2.getFitness() << "\nTO GET\n" 
			<< child.toString() << " " << child.getFitness() << "\n\n\n";
		std::cout << "========================================================================\n";
	}*/
	
	return child;
}

void GeneticBalancer::Chromosome::mutate()
{
	auto temp = *this;
	///randomly select eliminations
	std::vector<bool> willBeEliminated(genes.size());
	willBeEliminated[std::distance(genes.begin(),
									std::min_element(genes.begin(), genes.end(), [](std::vector<int> a, std::vector<int> b) { return a.size() < b.size(); }))
	                ] = true; //smallest bin
	if (parent.MAX_MUTATION_SEVERITY*genes.size() >= 3) //minimum 3 bins
	{
		for (int i = parent.random(2, parent.MAX_MUTATION_SEVERITY*genes.size() - 1); i > 0; --i)
		{
			int j;
			do { j = parent.random(0, genes.size() - 1); } while (willBeEliminated[j]);
			willBeEliminated[j] = true;
		}
	}
	//std::cout << "***********************for " << toString() << " will be eliminated ";
	//for (auto&x : willBeEliminated) std::cout << x;
	//std::cout << "\n";

	///eliminate
	std::vector<int> eliminated;
	for (int i = 0; i < willBeEliminated.size(); ++i)
	{
		if (willBeEliminated[i])
		{
			for (auto& item : genes[i]) eliminated.push_back(item);
			genes[i].clear();
		}
	}
	genes.erase(std::remove_if(
		genes.begin(), genes.end(),
		[](std::vector<int> x) {
		return x.empty();
	}), genes.end());
	
	///ff
	std::random_shuffle(eliminated.begin(), eliminated.end());
	parent.firstFit(eliminated, *this);

	///recalculate fitness
	calcFitness();
	/*
	if (getFitness() - temp.getFitness() > 0.1 * (1.0-temp.getFitness()))
	{
		std::cout << "*******************" << temp.toString() << temp.getFitness() << " mutated to \n";
		std::cout << "*******************" << toString() << getFitness() << "\n\n";
	}*/
}

void GeneticBalancer::Chromosome::inverse()
{
	//std::cout << "%%%%%%%%%%%%%%%%%%%%%Inversed\n" << toString() << " to\n";
	std::sort(genes.begin(), genes.end(),
		[](const std::vector<int>& a, const std::vector<int>& b) {
		return std::accumulate(a.begin(),a.end(),0) > std::accumulate(b.begin(), b.end(), 0);
	});
	//std::cout << toString() << "to\n";
	int n = genes.size();
	decltype(genes) temp(n);
	int tempI = 0;
	for (int i = n - 2; i >= 0; i -= 2) temp[tempI++] = genes[i];
	for (int i = !(n%2); i < n; i += 2) temp[tempI++] = genes[i];
	genes = temp;
	//std::cout << toString() << "\n";
}

std::string GeneticBalancer::Chromosome::toString() const
{
	std::stringstream ss;
	ss << "[";
	for (int i = 0; i < genes.size(); ++i)
	{
		for (int j = 0; j < genes[i].size(); ++j)
		{
			ss << genes[i][j];
			if (j < genes[i].size()-1) ss << ",";
		}
		int w = 23;
		if ((i + 1) * w + 1 > ss.str().size())
			ss << std::string((i + 1) * w + 1 - ss.str().size(), ' ');
		if (i < genes.size()-1) ss << "|";
	}
	ss << "]";
	return ss.str();
}

std::vector<std::vector<int>> GeneticBalancer::Chromosome::toBins()
{
	for (auto& workstation : genes)
		std::sort(workstation.begin(), workstation.end());
	std::sort(genes.begin(), genes.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });
	return genes;
}

//=============================================================================================================================================================
std::vector<std::vector<int>> GeneticBalancer::balance(std::vector<int> items, int binCapacity, const PrecedenceGraph& pg, std::vector<double>& bestFitness, long& elapsedTime)
{
	this->items = items;
	this->binCapacity = binCapacity;
	this->precedenceGraph = pg;
	srand(1337);

	std::clock_t start;
	start = std::clock();
	auto result = gga(bestFitness);
	elapsedTime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	return result;
}

/**
Genetic Grouping Algorithm

Population is always sorted by fitness
The first individual is the fittest
*/
std::vector<std::vector<int>> GeneticBalancer::gga(std::vector<double>& bestFitness)
{
	///init population
	auto population = initPopulation(POPULATION_SIZE);

	///evolution cycle
	for (int i = 0; i < MAX_NO_OF_GENERATIONS; ++i)
	{
		//printPopulation(population, i);
		//std::cout << "generation " << i << " best fitness: " << population[0].getFitness() << "\n";
		bestFitness.push_back(population[0].getFitness());
		///check population: if fitness is max -sucess, if diversity is minimal -failure
		if (population.front().isMaximallyFit() || population.front().getFitness() - population.back().getFitness() < DIVERSITY_THRESHOLD)
		{
			//std::cout << "\n[WARNING] BREAK ACCORDING TO DIVERSITY CRITERION\n";
			break;
		}

		///crossover
		decltype(population) children;
		for (int parentsCount = 0; parentsCount < int(POPULATION_SIZE*CROSSOVER_RATE) / 2; ++parentsCount)
		{
			auto parents = selectParents(population);
			children.push_back(crossover(parents.first, parents.second));
		}
		for (auto& individual : children)
			individual.inverse();
		population.insert(population.end(), children.begin(), children.end());

		///mutations
		decltype(population) mutants;
		for (auto& individual : population)
			if (randomZeroToOne() < MUTATION_RATE)
			{
				Chromosome mutant = individual;
				mutant.mutate();
				mutants.push_back(mutant);
			}
		for (auto& individual : mutants)
			individual.inverse();
		population.insert(population.end(), mutants.begin(), mutants.end());

		///prepare for the next generation
		sortPopulation(population);
		if (population.size() > POPULATION_SIZE)
			population.erase(population.begin()+POPULATION_SIZE, population.end());
	}
	//printPopulation(population, -1);
	return population.front().toBins();
}