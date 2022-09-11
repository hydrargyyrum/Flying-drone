#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>

sf::Color hsv2rgb(int hue)
{
    float sat = 1, val = 1.;
    hue %= 360;
    while (hue < 0) hue += 360;

    int h = hue / 60;
    float f = float(hue) / 60 - h;
    float p = val * (1.f - sat);
    float q = val * (1.f - sat * f);
    float t = val * (1.f - sat * (1 - f));

    switch (h)
    {
    default:
    case 0:
    case 6: return sf::Color(val * 255, t * 255, p * 255);
    case 1: return sf::Color(q * 255, val * 255, p * 255);
    case 2: return sf::Color(p * 255, val * 255, t * 255);
    case 3: return sf::Color(p * 255, q * 255, val * 255);
    case 4: return sf::Color(t * 255, p * 255, val * 255);
    case 5: return sf::Color(val * 255, p * 255, q * 255);
    }
}

class vertex : public sf::Vertex
{
public:
    sf::Vector3f pos;
    float scale = 6.;

    vertex() { }

    vertex(float x, float y, float z, sf::Color color) {
        this->pos = sf::Vector3f(x, y, z);
        this->color = color;
    }

    void rotate(double& angle)
    {
        double i, j, k, x;
        i = pos.x;
        j = pos.y;
        k = pos.z;
        x = i * std::cos(angle) - j * std::sin(angle);
        auto vect = sf::Vector2f(x, k);
        position = scale * vect + sf::Vector2f(400, 300);
    };

    vertex operator+(const vertex r) const
    {
        vertex l = *this;
        vertex res;
        res.color = l.color;
        res.pos = sf::Vector3f(
            l.pos.x + r.pos.x,
            l.pos.y + r.pos.y,
            l.pos.z + r.pos.z
        );
        return res;
    }

    void print() {
        std::cout << "x = " << pos.x << " y = " << pos.y << " z = " << pos.z << std::endl;
    }
};