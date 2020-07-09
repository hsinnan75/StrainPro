#include <unistd.h>
#include <sys/stat.h>
#include "structure.h"

int parseLine(char* line) {
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i - 3] = '\0';
	i = atoi(p);
	return i;
}

void ShowReadSequence(int64_t gPos)
{
	char *read = new char[ReadLength + 1];
	strncpy(read, RefSequence + gPos, ReadLength); read[ReadLength] = '\0';
	printf("%s\n", read); delete[] read;
}

string ShowTaxRank(int rank)
{
	switch (rank)
	{
	case 5: return "subspecies";
	case 10: return "species";
	case 15: return "subgenus";
	case 20: return "genus";
	case 25: return "subfamily";
	case 30: return "family";
	case 35: return "suborder";
	case 40: return "order";
	case 45: return "subclass";
	case 50: return "calss";
	case 55: return "subphylum";
	case 60: return "phylum";
	case 65: return "subkingdom";
	case 70: return "kingdom";
	default: return "root";
	}
}
