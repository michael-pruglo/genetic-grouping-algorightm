#include <SFML/Graphics.hpp>
#include "LineBalancingTester.h"
#include <iostream>

int main()
{
	sf::RenderWindow window(sf::VideoMode(1200, 650), "Line Balancer");
	window.setPosition(sf::Vector2i(1910-window.getSize().x, 1050 - window.getSize().y));

	int testNo = 1;
	window.clear();
	LineBalancingTester(&window).run(testNo);
	window.display();
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R)
			{
				window.clear();
				LineBalancingTester(&window).run(++testNo>5 ? (testNo=1) : testNo);
				window.display();
			}
		}
	}

	return 0;
}