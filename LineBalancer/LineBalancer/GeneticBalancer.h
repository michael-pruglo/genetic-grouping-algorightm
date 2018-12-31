#pragma once

#include <algorithm>
#include <vector>

class GeneticBalancer
{
	const int				POPULATION_SIZE = 10;
	const int				MAX_NO_OF_GENERATIONS = 3;
	const double			DIVERSITY_THRESHOLD = 0.01;
	const double			CROSSOVER_RATE = 0.2;
	const double			MUTATION_RATE = 0.1;
	const double			MAX_MUTATION_SEVERITY = 0.8;

	std::vector<int>		items;
	int						binCapacity;

public:
	class PrecedenceGraph
	{
		std::vector<std::vector<int>> longestPath;
	public:
		PrecedenceGraph() {};
		PrecedenceGraph(const std::vector<std::pair<int, int>>& edges);
		int getMaxDistance(int a, int b) { return longestPath[std::min(a, b)][std::max(a, b)]; };
	};
private:
	PrecedenceGraph						precedenceGraph;
	class Chromosome;

	Chromosome							randomChromosome();
	double								randomZeroToOne();
	int									random(int min, int max);
	int									spinRoulette(const std::vector<double>& probabilities);
	std::vector<Chromosome>				initPopulation(int size);
	std::pair<Chromosome, Chromosome>	selectParents(const std::vector<Chromosome>& population);
	std::vector<Chromosome>				selectParents2(const std::vector<Chromosome>& population);
	void								firstFit(const std::vector<int>& candidates, Chromosome& chromosome);
	Chromosome							crossover(const Chromosome& parent1, const Chromosome& parent2);
	void								sortPopulation(std::vector<Chromosome>& population);
	void								printPopulation(const std::vector<Chromosome>& population, int id);
	
	class Chromosome
	{
				GeneticBalancer&				parent;
				std::vector<std::vector<int>>	genes;
				double							fitness;

				void							calcFitness();

		friend	Chromosome						GeneticBalancer::randomChromosome();
	public:
				explicit Chromosome(GeneticBalancer& parent) : parent(parent) {}
				Chromosome& operator=(const Chromosome& b);
	public:
				double							getFitness() const { return fitness; }
				bool							isFitter(const Chromosome& b) const;
				bool							isMaximallyFit() const;
		friend	void							GeneticBalancer::firstFit(const std::vector<int>& candidates, Chromosome& chromosome);
		friend	Chromosome						GeneticBalancer::crossover(const Chromosome& parent1, const Chromosome& parent2);
				void							mutate();
				void							inverse();


				std::string						toString() const;
				std::vector<std::vector<int>>	toBins();
	};

private:
	std::vector<std::vector<int>> gga(std::vector<double>& bestFitness);
public:
	std::vector<std::vector<int>> balance(std::vector<int> items, int binCapacity, const PrecedenceGraph& pg, std::vector<double>& bestFitness, long& elapsedTime);
};

