#include <stdio.h>
#include <stdlib.h>

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ double enthält
double map[3][3];

// Die Funktion set_temperature soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
	void set_temperature (int x, int y, double temperature)
	{
		if(x >= 3 || x < 0 || y >= 3 || y < 0  )
		{
			fprintf(stderr, "x oder y müssen jeweils zwischen 0 und 3 liegen \n \n");
			return;
		}
		map[x][y] = temperature;
	}

	// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
	void show_map (void)
	{
		for(int x = 0; x < 3 ; x++)
		{
			for(int y = 0; y < 3; y++)
			{
				printf("%.2f\t",map[x][y]);
			}
			printf("\n");
		}
		printf("\n");
	}

// Die Funktion set_to_average soll an Position [x, y] den Durchschnitt der 8 umgebenen
// Temperaturen in das Array map eintragen.
// Für Werte außerhalb des Arrays soll der Wert 0 angenommen werden.
// Verwenden Sie hierfür auch die Funktion set_temperature.
void set_to_average (int x, int y)
{
    double res = 0.0;
    // row 1
    if( x-1 >= 0 && x-1 < 3 )
    {
        if( y-1 >= 0 && y-1 < 3 )
        {
            res += map[x-1][y-1];
        }
        if( y >= 0 && y < 3 )
        {
            res += map[x-1][y];
        }
        if( y+1 >= 0 && y+1 < 3 )
        {
            res += map[x-1][y+1];
        }
    }
    // row 2
    if( x >= 0 && x < 3 )
    {
        if( y-1 >= 0 && y-1 < 3 )
        {
            res += map[x][y-1];
        }
        if( y+1 >= 0 && y+1 < 3 )
        {
            res += map[x][y+1];
        }
    }
    // row 3
    if( x+1 >= 0 && x+1 < 3 )
    {
        if( y-1 >= 0 && y-1 < 3 )
        {
            res += map[x+1][y-1];
        }
        if( y >= 0 && y < 3 )
        {
            res += map[x+1][y];
        }
			if( y+1 >= 0 && y+1 < 3 )
			{
				res += map[x+1][y+1];
        }
    }
    set_temperature(x, y, res/ 8);
}

// In dieser Funktion darf nichts verändert werden!
int main (void)
{
	set_temperature(0, 1, 40);
	set_temperature(1, 0, 160);
	set_temperature(1, 4, 75);
	set_temperature(1, 2, 80);
	set_temperature(2, 1, 120);

	show_map();

	set_temperature(0, 0, 20.5);
	set_temperature(0, 2, 14.8);
	set_temperature(0, 2, 22.7);
	set_temperature(2, 0, 100.2);
	set_temperature(2, 2, 20.6);
	set_temperature(2, 2, 200.2);
	set_temperature(1, 3, 201.06);
	set_temperature(1, 1, 50.5);

	show_map();

  set_to_average(0,0);
  set_to_average(2,0);
  set_to_average(1,2);

  show_map();

	return 0;
}
