#include "structure.h"

#define BufferSize 1024

bool CheckBWAIndexFiles(char* IndexPrefix)
{
	fstream file;
	bool bChecked = true;

	file.open(((string)IndexPrefix + ".ann").c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	file.open(((string)IndexPrefix + ".amb").c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	file.open(((string)IndexPrefix + ".pac").c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	return bChecked;
}

bool CheckReadFormat(const char* filename)
{
	char buf[1];
	gzFile file = gzopen(filename, "rb");
	gzread(file, buf, 1); gzclose(file);

	if (buf[0] == '@') return true; // fastq
	else return false;
}

int IdentifyHeaderBoundary(char* str, int len)
{
	int i;

	for (i = 0; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '/' || str[i] == '\t') return i;
	}
	return len - 1; // remove the last character of '\n'
}

ReadItem_t GetNextEntry(FILE *file)
{
	ssize_t len;
	size_t size = 0;
	ReadItem_t read;
	char *buffer = NULL;

	read.seq = NULL; read.rlen = 0; read.AlnSummary.score = 0;

	if ((len = getline(&buffer, &size, file)) != -1)
	{
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, buffer + 1, len); read.header[len] = '\0';
		if (FastQFormat)
		{
			if ((read.rlen = getline(&buffer, &size, file)) != -1)
			{
				read.seq = new char[read.rlen];
				strncpy(read.seq, buffer, read.rlen);
				getline(&buffer, &size, file); getline(&buffer, &size, file);
				read.rlen -= 1; read.seq[read.rlen] = '\0';
			}
			else read.rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read.rlen = (int)seq.length()) > 0)
			{
				read.seq = new char[read.rlen + 1];
				strcpy(read.seq, (char*)seq.c_str());
				read.seq[read.rlen] = '\0';
			}
		}
	}
	if (buffer) free(buffer);

	return read;
}

int GetNextChunk(FILE *file, ReadItem_t* ReadArr)
{
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;
		if (++iCount == ReadChunkSize) break;
	}
	return iCount;
}

ReadItem_t gzGetNextEntry(gzFile file)
{
	ReadItem_t read;
	char buffer[BufferSize];

	read.seq = NULL; read.rlen = 0; read.AlnSummary.score = 0;

	if (gzgets(file, buffer, BufferSize) != NULL)
	{
		//len = IdentifyHeaderBoundary(buffer, strlen(buffer));
		if ((buffer[0] == '@' || buffer[0] == '>'))
		{
			//len -= 1; read.header = new char[len + 1];
			//strncpy(read.header, buffer + 1, len); read.header[len] = '\0';
			gzgets(file, buffer, BufferSize); read.rlen = strlen(buffer) - 1; read.seq = new char[read.rlen + 1]; read.seq[read.rlen] = '\0';
			strncpy(read.seq, buffer, read.rlen);

			if (FastQFormat)
			{
				gzgets(file, buffer, BufferSize);
				gzgets(file, buffer, BufferSize);
			}
		}
	}
	return read;
}

int gzGetNextChunk(gzFile file, ReadItem_t* ReadArr)
{
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		if (++iCount == ReadChunkSize) break;
	}
	return iCount;
}
