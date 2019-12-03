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

void NR_Analysis()
{
	FILE* fp;
	int64_t ori_size, nr_size;
	DIR *dir;
	int p;
	struct dirent *ds;
	string path, filename;

	dir = opendir(OutputFolder.c_str());

	ori_size = nr_size = 0;
	while ((ds = readdir(dir)) != NULL)
	{
		filename = ds->d_name; p = filename.find_last_of('.');
		if (filename.substr(p + 1) == "fna")
		{
			fp = fopen(ds->d_name, "r"); fseek(fp, 0, SEEK_END);
			ori_size += ftell(fp);
			fclose(fp);
		}
		else if (filename.substr(p + 1) == "nrs")
		{
			fp = fopen(ds->d_name, "r"); fseek(fp, 0, SEEK_END);
			nr_size += ftell(fp);
			fclose(fp);
		}
	}
	closedir(dir);

	printf("Ori_Size=%lld, NR_Size=%lld (%.4f)\n", (long long)ori_size, (long long)nr_size, (1.0*nr_size / ori_size));
}

void RemoveSeqFiles()
{
	int p;
	DIR *dir;
	char cmd[1024];
	struct dirent *ds;
	string filename, filetype;

	dir = opendir(OutputFolder.c_str());
	while ((ds = readdir(dir)) != NULL)
	{
		filename = ds->d_name; p = filename.find_last_of('.'); filetype = filename.substr(p + 1);
		if (filetype == "rep" || filetype == "clr1" || filetype == "clr2")
		{
			sprintf(cmd, "rm %s/%s", OutputFolder.c_str(), ds->d_name); system(cmd);
		}
	}
	closedir(dir);
}

void Remove_BWT_Files()
{
	int p;
	DIR *dir;
	char cmd[1024];
	struct dirent *ds;
	string filename, filetype;

	dir = opendir(OutputFolder.c_str());
	while ((ds = readdir(dir)) != NULL)
	{
		filename = ds->d_name; p = filename.find_last_of('.'); filetype = filename.substr(p + 1);
		if (filetype == "amb" || filetype == "ann" || filetype == "bwt" || filetype == "pac" || filetype == "sa")
		{
			sprintf(cmd, "rm %s/%s", OutputFolder.c_str(), ds->d_name); system(cmd);
		}
	}
	closedir(dir);
}
