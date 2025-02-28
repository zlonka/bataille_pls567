// original code : https://github.com/cristal-smac/bataille/blob/main/records_C4_naturelle.py
// records : https://github.com/cristal-smac/bataille/blob/main/records_C4_naturelle.py
// to compile on Linux : gcc -O4 bataille.cpp

#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define uchar unsigned char
#define uint unsigned int

#define OPT 1
#define NAT 2
#ifndef RANGEMENT
#define RANGEMENT OPT
#endif
#if RANGEMENT == OPT
#define RANGEMENT_STR "optimise"
#define FRECORDS "records.o.txt"
#define play playOpt
#elif RANGEMENT == NAT
#define RANGEMENT_STR "naturel"
#define FRECORDS "records.txt"
#define play playNat
#else
#error "RANGEMENT must be OPT or NAT."
#endif

#ifdef __unix__
#define RND lrand48()
#define SRAND srand48
#elif
#define RND rand()
#define SRAND srand
#endif

#define TYP uint

typedef struct{
    TYP cartes[52];
    int nb;
}JEU;
JEU j1, j2, centre;

typedef struct{
    JEU jeu;
    int score;
}PEOPLE;

char startFile[512];

uint record=0, loopMinRecord=9999;
void recordLoad(const int v) {
    FILE* f = fopen(FRECORDS, "r");
    if (f) {
        char buf[512];
        while (fgets(buf, 512, f)) {
            if (buf[0] == '#') continue;
            uint x, y;
            sscanf(buf, "%u", &x);
            if (x != v) continue;
            char kw[256];
            sscanf(buf, "%u %s", &x, kw);
            if (!strcmp(kw, "len")){
				sscanf(buf, "%u %s %u", &x, kw, &y);
				if (y > record) record = y;
			}
			else if (!strcmp(kw, "loop")){
				sscanf(buf, "%u %s %u", &x, kw, &y);
				if (y < loopMinRecord) loopMinRecord = y;
			}
        }
        fclose(f);
    }
}
void recordSave(const int v, const int score, const char *typ, const char *str) {
	if (!strcmp(typ, "len") || !strcmp(typ, "loop")){
		// printf("**************************************************\n");
		printf("***** NEW %s FOR v=%d : %d *****", typ, v, score);
		// printf("\n**************************************************\n");
	}
	else{
		fprintf(stderr, "unknown type '%s' to save, exit.\n", typ);		// ???
		exit(7);
	}
    // printf("saving %s record for v=%d : %d...", typ, v, score);
    printf("saving...");
    FILE* f = fopen(FRECORDS, "a");
    if (f) {
        fprintf(f, "%d\t%s\t%d\t%s", v, typ, score, str);
        time_t t = time(NULL);
        struct tm tim = *localtime(&t);
        fprintf(f, " %02d/%02d/%d %02d:%02d:%02d", tim.tm_mday, tim.tm_mon + 1, tim.tm_year + 1900, tim.tm_hour, tim.tm_min, tim.tm_sec);
        fprintf(f, "\n");
        fclose(f);
        printf("done.\n");
    }
    else printf(" ! ERROR WHILE SAVING.\n");
}

clock_t start;

void jeuDisp(JEU* j, int hilight=-1);
static int cmpCartes(const void *a, const void *b){
	if (*(TYP*)a < *(TYP*)b) return 1; else if (*(TYP*)a > *(TYP*)b) return -1; return 0;
}
__inline void jeuAddCarte(JEU* j, const TYP v) {
    j->cartes[j->nb++] = v;
}
void jeuAddCartes(JEU* j, const TYP* t, const int n) {
    for (int i = 0; i < n; i++) j->cartes[j->nb++] = t[i];
}
void jeuSetCartes(JEU* j, const TYP* t, const int n) {
    memcpy(j->cartes, t, n * sizeof(TYP)); j->nb = n; // for (int i = 0; i < n; i++) j->cartes[j->nb++] = t[i];
}
__inline void jeuClear(JEU* j) {
    j->nb = 0;
}
void jeuAddJeu(JEU* j, JEU* add) {
    for (int i = 0; i < add->nb; i++) j->cartes[j->nb++] = add->cartes[i];
    jeuClear(add);
}
// j+=add
__inline void jeuAddReverseJeu(JEU* j, JEU* add) {
    for (int i = 0; i < add->nb; i++) j->cartes[j->nb++] = add->cartes[add->nb - i - 1];
    jeuClear(add);
}
// sort(add), j+=add
__inline void jeuAddSortedJeu(JEU* j, JEU* add) {
	// printf("add="); jeuDisp(add); printf("\n");
    qsort(add->cartes, add->nb, sizeof(TYP), cmpCartes);
	// printf("add="); jeuDisp(add); printf("\n");
    // exit(1);
    for (int i = 0; i < add->nb; i++) j->cartes[j->nb++] = add->cartes[i];
    jeuClear(add);
}
TYP jeuPick(JEU* j) {
    TYP v = j->cartes[0];
    j->nb--;
    if (j->nb) {
        // for (int i = 0; i < j->nb; i++) j->cartes[i] = j->cartes[i + 1];
        memcpy(j->cartes, j->cartes + 1, j->nb * sizeof(TYP));
    }
    return v;
}
void sample1() {
    TYP t1[16] = { 3,2,1,1,2,1 };
    // jeuClear(&j1);    jeuAddCartes(&j1, t1, 6);
    jeuSetCartes(&j1, t1, 6);

    TYP t2[16] = { 2,1,3,2,3,3 };
    // jeuClear(&j2);    jeuAddCartes(&j2, t2, 6);
    jeuSetCartes(&j2, t2, 6);

    jeuClear(&centre);
}
void sample2() {
    TYP t1[16] = { 3,1,1,2,3,2 };
    // jeuClear(&j1);    jeuAddCartes(&j1, t1, 6);
    jeuSetCartes(&j1, t1, 6);

    TYP t2[16] = { 3,1,3,2,1,2 };
    // jeuClear(&j2);    jeuAddCartes(&j2, t2, 6);
    jeuSetCartes(&j2, t2, 6);

    jeuClear(&centre);
}
void distrib(JEU *j) {
    int n = j->nb >> 1;
    /*
    jeuClear(&j1);    jeuAddCartes(&j1, j->cartes, n);
    jeuClear(&j2);    jeuAddCartes(&j2, j->cartes + n, n);
    */
    jeuSetCartes(&j1, j->cartes, n);
    jeuSetCartes(&j2, j->cartes + n, n);
    jeuClear(&centre);
}
void jeuDisp(JEU* j, int hilight) {
    printf("[");
    for (int i = 0; i < j->nb; i++) {
        if (i > 0) printf(",");
        if (i == hilight) printf("*%d*", j->cartes[i]); else printf("%d", j->cartes[i]);
    }
    printf("]");
}
void jeuDispFile(FILE *f, JEU* j, int hilight=-1) {
    fprintf(f, "[");
    for (int i = 0; i < j->nb; i++) {
        if (i > 0) fprintf(f, ",");
        if (i == hilight) fprintf(f, "*%d*", j->cartes[i]); else fprintf(f, "%d", j->cartes[i]);
    }
    fprintf(f, "]");
}
void jeuDispStr(JEU* j, char *str) {
    strcpy(str, "[");
    for (int i = 0; i < j->nb; i++) {
        if (i > 0) strcat(str, ",");
        char word[80];
        sprintf(word, "%d", j->cartes[i]);
        strcat(str, word);
    }
    strcat(str, "]");
}
void jeuCopyTo(JEU *jeu, JEU *jeuCopy) {
    jeuCopy->nb = jeu->nb;
    memcpy(jeuCopy->cartes, jeu->cartes, jeu->nb * sizeof(TYP));
}
void jeuShuffle(JEU* j, const int nbPermut=0) {
    int nbP = nbPermut != 0 ? nbPermut : j->nb;
    for (int i = 0; i < nbP; i++) {
        int a = RND % j->nb, b = RND % j->nb;
        TYP tmp = j->cartes[a]; j->cartes[a] = j->cartes[b]; j->cartes[b] = tmp;
    }
}
void disp() {
    /*
    printf("J1 : ");
    jeuDisp(&j1);
    printf(" J2 : ");
    jeuDisp(&j2);
    printf(" Centre : ");
    jeuDisp(&centre);
    */
    // printf("\n");
    printf("J1/J2/C : ");
    for (int i = 0; i < j1.nb; i++) {
        if (i > 0) printf(",");
        printf("%d", j1.cartes[i]);
    }
    printf(" / ");
    for (int i = 0; i < j2.nb; i++) {
        if (i > 0) printf(",");
        printf("%d", j2.cartes[i]);
    }
    printf(" / ");
    for (int i = 0; i < centre.nb; i++) {
        if (i > 0) printf(",");
        printf("%d", centre.cartes[i]);
    }
}

int iRound = 0, roundLimit = 3000, globNbCartes = 0;

TYP playNat(bool trace = true) {
    for (iRound = 0;; iRound++) {
        if (iRound > roundLimit) return 3;
        if (trace) { printf("Pli %04d ", iRound);  disp(); }
        if (j1.nb == 0 && j2.nb == 0) return 0; // deuce
        if (j1.nb == 0) return 2;   // 2 wins
        if (j2.nb == 0) return 1;   // 1 wins
        TYP c1 = jeuPick(&j1), c2 = jeuPick(&j2);
        if (c1 > c2) {
            if (trace) printf("J1 wins ");
            jeuAddCarte(&j1, c1);
            jeuAddCarte(&j1, c2);
            jeuAddReverseJeu(&j1, &centre);
        }
        else if (c2 > c1) {
            if (trace) printf("J2 wins ");
            jeuAddCarte(&j2, c2);
            jeuAddCarte(&j2, c1);
            jeuAddReverseJeu(&j2, &centre);
        }
        else {
            if (trace) printf("Bataille ");
            jeuAddCarte(&centre, c1);
            jeuAddCarte(&centre, c2);
            // if (trace) { jeuDisp(&centre); }
        }
        if (trace) { printf(" => ");  disp(); printf("\n"); }
    }
}
typedef struct{
	JEU j1, j2, c;
}LOOP;
LOOP loop[8192];
int detectLoop(const int nb){
	for(int i=nb-1;i>=0;i--){
		if(loop[i].j1.nb==j1.nb && loop[i].j2.nb==j2.nb && loop[i].c.nb==centre.nb
			&& !memcmp(loop[i].j1.cartes, j1.cartes, j1.nb * sizeof(TYP))
			&& !memcmp(loop[i].j2.cartes, j2.cartes, j2.nb * sizeof(TYP))
			&& !memcmp(loop[i].c.cartes, centre.cartes, j1.nb * sizeof(TYP))
		){
			return nb-i;
		}
	}
	return 0;
}
int iLoopMax=0, loopSizeMin=999999;
void addLoop(const int i){
	if (i > iLoopMax){
		iLoopMax = i;
		// fprintf(stderr, "iLoopMax=%d\n", iLoopMax);
	}
	memcpy(&loop[i].j1, &j1, sizeof(JEU));
	memcpy(&loop[i].j2, &j2, sizeof(JEU));
	memcpy(&loop[i].c, &centre, sizeof(JEU));
}
TYP playOpt(bool trace = true) {
    iRound = 0;
    for (int iRoundPli = 0;; iRound++) {
        if (iRound > roundLimit) return 3;

        int loopSize = detectLoop(iRoundPli-1);
        if (loopSize > 0){
			/*
			if (loopSize < loopSizeMin){
				// loopSizeMin = loopSize;
				printf("new min loop size : %d\n", loopSize);
				printf("iRound=%d loop size %d detected\n", iRound, loopSize);
				for(int i=0;i<iRound;i++){ printf("%4d", i); jeuDisp(&loop[i].j1); printf(" "); jeuDisp(&loop[i].j2); printf(" "); jeuDisp(&loop[i].c); printf("\n"); }
				int k = iRound-1; printf("%4d", k); jeuDisp(&loop[k].j1); printf(" "); jeuDisp(&loop[k].j2); printf(" "); jeuDisp(&loop[k].c); printf("\n");
				k -= loopSize; printf("%4d", k); jeuDisp(&loop[k].j1); printf(" "); jeuDisp(&loop[k].j2); printf(" "); jeuDisp(&loop[k].c); printf("\n");
				// exit(1);
			}
			*/
			iRound = loopSize;
			return 3;
        }

        if (trace) { printf("Pli %04d ", iRound);  disp(); }
        if (j1.nb == 0 && j2.nb == 0) return 0; // deuce
        if (j1.nb == 0) return 2;   // 2 wins
        if (j2.nb == 0) return 1;   // 1 wins
        TYP c1 = jeuPick(&j1), c2 = jeuPick(&j2);
        if (c1 > c2) {
            if (trace) printf("J1 wins ");
            // jeuAddCarte(&j1, c1); jeuAddCarte(&j1, c2); jeuAddReverseJeu(&j1, &centre);
            jeuAddCarte(&centre, c2);
            jeuAddCarte(&centre, c1);
            jeuAddSortedJeu(&j1, &centre);
	        addLoop(iRoundPli++);
        }
        else if (c2 > c1) {
            if (trace) printf("J2 wins ");
            // jeuAddCarte(&j2, c2); jeuAddCarte(&j2, c1); jeuAddReverseJeu(&j2, &centre);
            jeuAddCarte(&centre, c1);
            jeuAddCarte(&centre, c2);
            jeuAddSortedJeu(&j2, &centre);
	        addLoop(iRoundPli++);
        }
        else {
            if (trace) printf("Bataille ");
            jeuAddCarte(&centre, c1);
            jeuAddCarte(&centre, c2);
            // if (trace) { jeuDisp(&centre); }
        }
        if (trace) { printf(" => ");  disp(); printf("\n"); }
    }
}
int cleanupCfg(const char* cfg, JEU *dst) {
    dst->nb = 0;
    int chk4[256];
    for (uint i = 1; i < 256; i++) chk4[i]=0;
    for (uint i = 0; i < strlen(cfg);) {
        TYP v = 0;
        while (cfg[i] >= '0' && cfg[i] <= '9') {
            v *= 10; v += cfg[i] - '0';
            i++;
        }
        if (v != 0) {
            dst->cartes[dst->nb] = v;
            if (v<1 || v>=256){
				fprintf(stderr, "Error : v=%d is not a correct value for a card (not in [1..256])\n", v);
				exit(5);
			}
			chk4[v]++;
            dst->nb++;
        }
        while (i < strlen(cfg) && (cfg[i] > '9' || cfg[i] < '0')) {
            // printf("ignore char '%c'\n", cfg[i]);
            i++;
        }
    }
    for (uint i = 1; i <= dst->nb / 4; i++){
		if (chk4[i]!=4){
			fprintf(stderr, "Error : %d has %d cards instead of 4\n", i, chk4[i]);
			exit(5);
		}
	}
    return dst->nb;
}
bool jeuOK(JEU* j, const TYP vmax) {  // check 4 cards of each v=1..vmax
    int occ[16];
    for (TYP v = 1; v <= vmax; v++) occ[v] = 0;
    for (int i = 0; i < j->nb;i++) {
        occ[j->cartes[i]]++;
        if (occ[j->cartes[i]] > 4) {
            // if (j->cartes[i] > 2) { jeuDisp(j); printf(" nb(%d) > 4\n", j->cartes[i]); }
            return false;
        }
    }
    // for (TYP v = 1; v <= vmax; v++) if (occ[v] != 4) return false;
    return true;
}
int roundMax = 0;
unsigned long long nbCfg = 0LL;
int nbv[16];
bool dispInfinite = false;
void recurInit(const TYP vmax) {
    for (TYP i = 1; i <= vmax; i++) nbv[i] = 0;
}
void recur(const int prof, const TYP vmax, JEU* j) {
    if (prof >= j->nb) {
        if (true || jeuOK(j, vmax)) {
            // jeuDisp(j);
            nbCfg++;
            distrib(j);
            TYP win = play(false);
            if (win == 0) {
                // printf("Deuce\n");
            }
            else if (win == 3) {
                if (dispInfinite) {
                    distrib(j);
                    printf("J1 : ");
                    jeuDisp(&j1);
                    printf(" J2 : ");
                    jeuDisp(&j2);
                    printf("Infinite loop (>= %d loops)\n", roundLimit);
                }
            }
            else {
                // printf("*** winner : %d iRound=%d\n", win, iRound);
                if (iRound > roundMax) {
                    roundMax = iRound;
                    // jeuDisp(&jeu); printf(" :: ");
                    distrib(j);
                    printf("J1 : ");
                    jeuDisp(&j1);
                    printf(" J2 : ");
                    jeuDisp(&j2);
                    clock_t stop = clock();
                    printf(" : new nbRoundMax=%d/%d (%.3f s.)\n", roundMax, roundLimit, (double)(stop - start) / CLOCKS_PER_SEC);
                }
            }
        }
        return;
    }
    for (TYP v = 1; v <= vmax; v++) {
        if (nbv[v] == 4) continue;
        nbv[v]++;
        j->cartes[prof] = v;
        recur(prof + 1, vmax, j);
        nbv[v]--;
    }
}
unsigned long long nbCfgDepart = 0LL;
void endGame(JEU* j, const TYP vmax) {
    nbCfg++;
    if (nbCfg < nbCfgDepart) return;
    if (!(nbCfg & 0xffff)){	// only update every 65536
		printf("cfg #%llu ", nbCfg);
		jeuDisp(j);
		printf(" saving...");
		FILE *f=fopen(startFile, "w");
		if (f==(FILE*)NULL){
			printf("error saving %s", startFile);
		}
		else{
			fprintf(f, "%llu", nbCfg);
			fclose(f);
			printf("done.");
		}
		printf("   %c", 13);
	}
    distrib(j);
    TYP win = play(false);
    if (win == 0) {
        // printf("Deuce\n");
    }
    else if (win == 3) {
        if (dispInfinite) {
            distrib(j);
            printf("J1 : ");
            jeuDisp(&j1);
            printf(" J2 : ");
            jeuDisp(&j2);
            printf("Infinite loop (>=%d)\n", roundLimit);
        }
    }
    else {
        // printf("*** winner : %d iRound=%d\n", win, iRound);
        if (iRound > roundMax) {
            roundMax = iRound;
            // jeuDisp(&jeu); printf(" :: ");
            distrib(j);
            printf("J1 : ");
            jeuDisp(&j1);
            printf(" J2 : ");
            jeuDisp(&j2);
            clock_t stop = clock();
            printf(" : new nbRoundMax=%d/%d loop #%llu (%.3f s.)\n", roundMax, roundLimit, nbCfg, (double)(stop - start) / CLOCKS_PER_SEC);
            if (roundMax > record) {
                char buf[512], word[256];
                jeuDispStr(&j1, buf);
                strcat(buf, " ");
                jeuDispStr(&j2, word);
                strcat(buf, word);
                recordSave((uint)vmax, roundMax, "len", buf);
                roundLimit = roundMax + 100;
            }
        }
    }
}

void allRecur(JEU* jeu, const TYP v, const TYP vmax) {
	// place 4 cartes v
    const int nb = jeu->nb;
    // if (jeu->cartes[0] > vmax/2) return;
    if (v == vmax) {
        endGame(jeu, vmax);
        return;
    }
    for (int i = 0; i < nb - 3; i++) {
        if (jeu->cartes[i] != vmax) continue;
        jeu->cartes[i] = v;
        for (int j = i + 1; j < nb - 2; j++) {
            if (jeu->cartes[j] != vmax) continue;
            jeu->cartes[j] = v;
            for (int k = j + 1; k < nb - 1; k++) {
                if (jeu->cartes[k] != vmax) continue;
                jeu->cartes[k] = v;
                for (int l = k + 1; l < nb; l++) {
                    if (jeu->cartes[l] != vmax) continue;
                    jeu->cartes[l] = v;
                    allRecur(jeu, v + 1, vmax);
                    jeu->cartes[l] = vmax;
                }
                jeu->cartes[k] = vmax;
            }
            jeu->cartes[j] = vmax;
        }
        jeu->cartes[i] = vmax;
    }
}
int calcNbBitsLeft(const uint v, const int nSur2){
	int nbOnes = 0;
	uint msk = 1 << nSur2;
	for(int i=0;i<nSur2;i++){
		if (v & msk) nbOnes++;
		msk <<= 1;
	}
	return nbOnes;
}
int calcNbBitsRight(const uint v, const int nSur2){
	int nbOnes = 0;
	uint msk = 1;
	for(int i=0;i<nSur2;i++){
		if (v & msk) nbOnes++;
		msk <<= 1;
	}
	return nbOnes;
}
void allOne(JEU* jeu, const TYP vmax) {
    const int nb = jeu->nb;
    for (int i = 0; i < nb; i++) jeu->cartes[i] = vmax;
    int nbSur2 = nb >> 1;
    // place 4 cartes 1, les 2 premières sont forcement dans la 1ere moitie, sinon on parcourt des solutions symetriques
    //                   ex : xxxx11x1 = 1x11xxxx
    for (int i = 0; i < nbSur2 - 1; i++) {
        jeu->cartes[i] = 1;
        for (int j = i + 1; j < nbSur2; j++) {
            jeu->cartes[j] = 1;
            for (int k = j + 1; k < nb - 1; k++) {
                jeu->cartes[k] = 1;
                for (int l = k + 1; l < nb; l++) {
                    jeu->cartes[l] = 1;
					/*
                    uint msk = (1<<(nb-i-1)) | (1<<(nb-j-1)) | (1<<(nb-k-1)) | (1<<(nb-l-1));
                    int nbBitsLeft=calcNbBitsLeft(msk, nbSur2);
                    int nbBitsRight=calcNbBitsRight(msk, nbSur2);
                    if (nbBitsLeft < nbBitsRight) continue;
					*/
                    allRecur(jeu, 2, vmax);
                    jeu->cartes[l] = vmax;
                }
                jeu->cartes[k] = vmax;
            }
            jeu->cartes[j] = vmax;
        }
        jeu->cartes[i] = vmax;
    }
}
int tDep[4], nbCfgOne=0;
void singleOne(JEU* jeu, const TYP vmax) {
    const int nb = jeu->nb;
    for (int i = 0; i < nb; i++) jeu->cartes[i] = vmax;
    // place 4 cartes 1
    int i = tDep[0], j=tDep[1], k=tDep[2], l=tDep[3];
	jeu->cartes[i] = 1;
    jeu->cartes[j] = 1;
	jeu->cartes[k] = 1;
	jeu->cartes[l] = 1;
	printf("ONEs : %d %d %d %d nb=%d\n", i, j, k, l, nb);
	allRecur(jeu, 2, vmax);
}
void allOneV2(JEU* jeu, const TYP vmax) {
    const int nb = jeu->nb;
    for (int i = 0; i < nb; i++) jeu->cartes[i] = vmax;
    // place 4 cartes 1
    int i = tDep[0], j=tDep[1], k=tDep[2], l=tDep[3];
	jeu->cartes[i] = 1;
    jeu->cartes[j] = 1;
	jeu->cartes[k] = 1;
	jeu->cartes[l] = 1;
	for(;;){
		nbCfgOne++;
		printf("ONEs : %d %d %d %d nb=%d\n", i, j, k, l, nb);
		allRecur(jeu, 2, vmax);
		// next config
		l++;
		if (l == nb){
			k++;
			if (k == nb - 1){
				j++;
				if (j == nb - 2){
					i++;
					if (i == nb - 3) return;
					j = i + 1;
				}
				k = j + 1;
			}
			l = k + 1;
		}
	}
}
void newBestScore(const int v, JEU *j, const int iRound, const uint attempt, char *str=NULL){
	roundMax = iRound;
	distrib(j);
	printf("gen#%u :: ", attempt);
	printf("J1 : ");
	jeuDisp(&j1);
	printf(" J2 : ");
	jeuDisp(&j2);
	printf(" : new nbRoundMax=%d/%d", roundMax, roundLimit);
	if (str != (char*)NULL) printf(" %s", str);
	printf("\n");
	if (roundMax > record) {
		roundLimit = roundMax + 100;
		record = roundMax;
		char buf[512], word[256];
		jeuDispStr(&j1, buf);
		strcat(buf, " ");
		jeuDispStr(&j2, word);
		strcat(buf, word);
		sprintf(word, " after %u attempts", attempt);
		strcat(buf, word);
		if (str != (char*)NULL){
			strcat(buf, " ");
			strcat(buf, str);
		}
		recordSave((uint)v, roundMax, "len", buf);
	}
}
static int cmpPeople(const void *p1, const void *p2){
	int a=((PEOPLE*)p1)->score, b=((PEOPLE*)p2)->score;
	return a==b ? 0 : (a > b ? -1 : 1);
}
int r1n(const int n){	// return x in [1..n]
	if (n <= 1) return 1;
	return 1+(RND % n);
}
int game(JEU jeu, const uint attempt, const int v){
	distrib(&jeu);					// updates j1 and j2
	TYP win = play(false);		// plays j1 and j2, updates iRound
	if (win == 0) { // printf("Deuce\n");
		return 0;
	}
	else if (win == 3) { // printf("Inifinite loop\n");
		if (jeu.nb>=8*4){
			fprintf(stderr, "Loop detected !!!\n");
			char line[1024], word[512];
			strcpy(line, "Loop detected for game: ");
			jeuDispStr(&jeu, word);	strcat(line, word);
			FILE *f=fopen(FRECORDS, "a");
			if (f){
				fprintf(f, "# %s\n", line);
				fclose(f);
			}
			else{
				fprintf(stderr, "%s\n", line);
				fprintf(stderr, "(err: can't update file results.txt)\n");
			}
			exit(1);
		}
		return 0;
	}
	else {
		// printf("*** winner : %d iRound=%d\n", win, iRound);
		if (iRound > roundMax) {
			roundMax = iRound;
			// jeuDisp(&jeu); printf(" :: ");
			distrib(&jeu);
			printf("J1 : ");        jeuDisp(&j1);
			printf(" J2 : ");       jeuDisp(&j2);
			printf(" : new nbRoundMax=%d/%d", roundMax, roundLimit);
			printf(" after %u attempts", attempt);
			printf("\n");
			if (roundMax > record) {
				roundLimit = roundMax + 100;
				record = roundMax;
				char buf[512], word[256];
				jeuDispStr(&j1, buf);
				strcat(buf, " ");
				jeuDispStr(&j2, word);
				strcat(buf, word);
				sprintf(word, " after %u attempts", attempt);
				strcat(buf, word);
				recordSave((uint)v, roundMax, "len", buf);
			}
		}
		return iRound;
	}
}
void usage(char *progname, const char *str=NULL){
	fprintf(stderr, "usage: %s [-h] [v=<v>] [cfg=<config>] [calcOne] [startOne=<a,b,c,d>] [random=<method>] [seed=<seed>] [genetic]\n", progname);
	fprintf(stderr, " <config> = sample1 | sample2 | free string like '1,2,2,1,2,2,1,1' or '1 3 2 4 1 3 2 4 1 2 3 4 1 2 3 4'\n");
	fprintf(stderr, "  sample1 : j1={3,2,1,1,2,1} j2={2,1,3,2,3,3}\n");
	fprintf(stderr, "  sample2 : j1={3,1,1,2,3,2} j2={3,1,3,2,1,2}\n");
	fprintf(stderr, " <startOne> = index of 1s in the game, ex: 0,1,2,3 or 2,6,8,9\n");
	fprintf(stderr, " <method> = brute | montecarlo | mc2\n");
	fprintf(stderr, " <seed> = <integer> | time (default=time)\n");
	fprintf(stderr, " parentsPercentage=<int> [default=25]\n");
	fprintf(stderr, " mutationRate=<int> [default=10]");
	if (str) fprintf(stderr, "%s\n", str);
	exit(1);
}
void dispDelai(FILE *dst, time_t t0, time_t t1){
	time_t dt = t1-t0;
	int h = (int)dt/3600;
	int m = ( dt - h * 3600 ) / 60;
	int s = dt % 60;
	if (h != 0) fprintf(dst, "%dh%02d'%02d", h, m, s);
	else fprintf(dst, "%d'%02d", m, s);
}
uint mirrorBits(const uint v, const int nbBits){
	uint ret = 0, mskRead = 1, mskWrite = 1;
	mskRead <<= nbBits - 1;
	for(int i=0;i<nbBits;i++){
		if (v & mskRead) ret |= mskWrite;
		mskWrite <<= 1;
		mskRead >>= 1;
	}
	return ret;
}
void swapCartes(JEU *j, const int a, const int b){
	TYP tmp=j->cartes[a];
	j->cartes[a] = j->cartes[b];
	j->cartes[b] = tmp;
}
int tOcc[64];
void rFill(JEU* j, const int v, const int prof, const int profmax){
	if (prof >= profmax){
		// printf("rFill prof=%d profmax=%d ", prof, profmax); jeuDisp(j); printf("\n");
		endGame(j, v);
		return;
	}
	for(int i=1;i<=v;i++){
		if (tOcc[i]==4) continue;
		tOcc[i]++;
		j->cartes[prof]=i;
		// printf(" rFill : j->cartes[%d]=%d (tOcc[%d]=%d)\n", prof, i, i, tOcc[i]);
		rFill(j, v, prof+1, profmax);
		tOcc[i]--;
	}
}
int jeuFromStr(JEU *j, const int v, char *str){
	// converts "1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,9" to [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,9]
	char *ptr = str;
	int tPartial[128];
	int iv=0, vCur=0;
	while(*ptr){
		while(*ptr>='0' && *ptr<='9'){ vCur*=10;vCur+=(*ptr - '0'); ptr++; }
		while(*ptr && (*ptr < '0' || *ptr > '9')){ ptr++; }
		if (vCur == 0 || vCur > v){
			fprintf(stderr, "partial must be in [1..%d]. %d found.\n", v, vCur);
			exit(1);
		}
		if (iv >= j->nb || iv >= 128){
			fprintf(stderr, "partial string has too many values (%d values don't fit in min(%d,128))\n", iv, j->nb);
			exit(1);
		}
		tPartial[iv++]=vCur;
		vCur=0;
	}
	printf("tPartial=["); for(int i=0;i<iv;i++) printf("%d ", tPartial[i]); printf("]\n");
	if (iv==0){
		fprintf(stderr, "partial string is empty\n");
		exit(1);
	}
	for(int i=0;i<j->nb;i++) j->cartes[i] = 0;
	for(int i=0;i<iv;i++) j->cartes[i] = tPartial[i];
	return iv;
}
void affDelai(long int dt, FILE *f=stdout){
	bool disp=false;
	if (dt>=3600*24){ long int n=dt/(3600*24); fprintf(f, "%ldd", n); dt-=n*(3600*24); disp=true; }
	if (dt>=3600 || disp){ long int n=dt/3600; fprintf(f, "%ldh", n); dt-=n*3600; disp=true; }
	if (dt>=60 || disp){ long int n=dt/60; if (disp) fprintf(f, "%02ldm", n); else fprintf(f, "%ldm", n); dt-=n*60; disp=true; }
	if (disp) fprintf(f, "%02ld'", dt); else fprintf(f, "%ld'", dt);
}
int main(int argc, char** argv) {
    int c = 4;  // nombre de couleurs
    int v = 3;  // nombre de valeurs
    char cfg[512] = "nope";
    char randomMethod[512] = "---";
    char partial[512];
    int vDep = 8;
	sprintf(startFile, "start.%d.txt", vDep);
    bool randomOn = false, geneticOn = false, singleOneOn = false, calcOneOn = false, partialOn = false;
	uint seed = 0;
	int parentsPercentage = 25;
	int mutationRatePer1000 = 100;

	// not to start from scratch
	tDep[0] = 0; tDep[1] = 1; tDep[2] = 2;  tDep[3] = 3;
	// default values
	strcpy(partial, "");

    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (!strncmp(argv[i], "v=", 2)) {
                sscanf(argv[i] + 2, "%d", &vDep);
				sprintf(startFile, "start.%d.txt", vDep);
            }
            else if (!strncmp(argv[i], "cfg=", 4)) {
                strcpy(cfg, argv[i] + 4);
            }
            else if (!strncmp(argv[i], "calcOne", 7)) {
                calcOneOn = true;
            }
            else if (!strncmp(argv[i], "partial=", 8)) {
                partialOn = true;
                strcpy(partial, argv[i]+8);
            }
            else if (!strncmp(argv[i], "depart=", 7)) {
                strcpy(partial, argv[i]+7);
            }
            else if (!strncmp(argv[i], "startOne=", 9)) {
                sscanf(argv[i] + 9, "%d,%d,%d,%d", tDep, tDep+1, tDep+2, tDep+3);
                singleOneOn = true;
            }
            else if (!strncmp(argv[i], "random=", 7)) {
                randomOn = true;
                strcpy(randomMethod, argv[i] + 7);
            }
            else if (!strncmp(argv[i], "seed=", 5)) {
                if (strcmp(argv[i] + 5, "time")) sscanf(argv[i] + 5, "%u", &seed);
            }
            else if (!strncmp(argv[i], "genetic", 7)) {
                geneticOn = true;
            }
            else if (!strncmp(argv[i], "parentsPercentage=", 18)) {
                sscanf(argv[i] + 18, "%d", &parentsPercentage);
            }
            else if (!strncmp(argv[i], "mutationRate=", 13)) {
                int tmp;
                sscanf(argv[i] + 13, "%d", &tmp);
                mutationRatePer1000 = tmp * 10;
            }
            else if (!strncmp(argv[i], "start=", 6)) {
                if (!strncmp(argv[i]+6, "auto", 4)){
					FILE *f=fopen(startFile, "r");
					if (f==(FILE*)NULL){
						fprintf(stderr, "Ooops, start=auto, but can't read file %s\n", startFile);
						exit(2);
					}
					int nbScan = fscanf(f, "%llu", &nbCfgDepart);
					fclose(f);
				}
				else{
					sscanf(argv[i] + 6, "%llu", &nbCfgDepart);
				}
				printf("starting at %llu\n", nbCfgDepart);
            }
            else if (!strncmp(argv[i], "-h", 2)) {
				usage(argv[0]);
            }
            else {
                fprintf(stderr, "option %s unknown\n", argv[i]);
                usage(argv[0]);
            }
        }
        if (randomOn){
			if (strcmp(randomMethod, "montecarlo")
				&& strcmp(randomMethod, "mc2")
				&& strcmp(randomMethod, "brute")
				)
			usage(argv[0], " random=[brute|montecarlo|mc2]");
		}
    }

    if (vDep > 13) {
        fprintf(stderr, "Sorry, v=%d too big, v must be in [1,13].\n", vDep);
        exit(1);
    }

    v = vDep;
    int nbCartes = v * c;

	if (seed==0) seed = time(NULL);
    SRAND(seed);

	JEU jDep;
    bool doSample = false;
    if (strcmp(cfg, "nope") /*!strcmp(cfg, "sample1") || !strcmp(cfg, "sample2")*/) doSample = true;

    TYP valeur = 1;
    v = vDep;
    if (v >= 4) record = 190;
    if (v >= 5) record = 562;
    if (v >= 6) record = 178;
    if (v >= 7) record = 1768;
    if (v >= 8) record = 2546;
    if (v >= 9) record = 2926;
    if (v >= 10) record = 3606;
    if (v >= 11) record = 1200;
    if (v >= 12) record = 4860;
    if (v >= 13) record = 5942;
    recordLoad(v);
    roundLimit = record + 100;
    printf("*** MODE RANGEMENT : %s *** current record : %u, loopMinRecord : %u, roundLimit=%d\n", RANGEMENT_STR, record, loopMinRecord, roundLimit);
    globNbCartes = nbCartes = v * c;

    if (doSample) {
        // distrib
        if (!strcmp(cfg, "sample1")) sample1();
        else if (!strcmp(cfg, "sample2")) sample2();
        else{
			JEU jeu;
			int nbCfg = cleanupCfg(cfg, &jeu);
            if (nbCfg != nbCartes) {
                fprintf(stderr, "bad cfg length : %d insteads of %d\n", nbCfg, nbCartes);
                return 1;
            }
            distrib(&jeu);
		}
        // disp
        disp();
        // play;
        printf("*** Playing !\n");
        TYP win = play(true);
        if (win == 0) printf("Deuce");
        else if (win == 3) printf("Inifinite loop (>=%d)", roundLimit);
        else printf("*** winner : %d", win);
        printf("\n");
        exit(0);
    }

    JEU jeu;
    // JEU jeuCopy;
    // init jeu
    jeu.nb = nbCartes;

    if (calcOneOn){
		bool disp = true;
		uint nbCases = 0, nbCasesTot = 0;
		int nbSur2 = nbCartes/2;
		uint mskLow = 0;
		char fOne[80];
		sprintf(fOne, "one.%d.txt", v);
		FILE *f = fopen(fOne, "w");
		if (f == (FILE*)NULL){
			fprintf(stderr, "Ooops, can't create file %s.\n", fOne);
			exit(2);
		}
		for(int i=0;i<nbSur2;i++) mskLow |= (1 << i);
		printf("v=%d nbCartes=%d nbSur2=%d mskLow=0x%08x\n", v, nbCartes, nbSur2, mskLow);
		for(int i=0;i<nbSur2;i++){
			for(int j=i+1;j<nbSur2;j++){
				uint vUp = (1<<(nbSur2-i-1)) | (1<<(nbSur2-j-1));
				if (disp) printf("vUp=0x%02x \n", vUp);
				for(int k=j+1;k<nbCartes-1;k++){
					for(int l=k+1;l<nbCartes;l++){
						nbCasesTot++;
						if (k>=nbSur2){
							uint vDown = (1<<(nbCartes-k-1)) | (1<<(nbCartes-l-1));
							if (disp) printf(" vDown= 1 << %d | 1 << %d = 0x%04x ", (nbCartes-k-1), (nbCartes-l-1), vDown);
							vDown &= mskLow;
							vDown = mirrorBits(vDown, nbSur2);
							if (disp) printf("i=%d j=%d k=%d l=%d vUp=0x%04x vDown=0x%04x ", i, j, k, l, vUp, vDown);
							if (vDown > vUp){
								if (disp) printf(" > REJECT\n");
								continue;
							}
						}
						if (disp) printf("i=%d j=%d k=%d l=%d > OK\n", i, j, k, l);
						fprintf(f, "%d,%d,%d,%d\n", i, j, k, l);
						// if (k==nbSur2 && l==32) exit(1);
						nbCases++;
					}
				}
			}
		}
		fclose(f);
        fprintf(stderr, "Finished. %u/%u cases. File %s created.\n", nbCases, nbCasesTot, fOne);
        exit(0);
	}
	if (singleOneOn){
        printf("singleOneOn (starting with [");
        for(int i=0;i<4;i++) if (i!=0) printf(",%d", tDep[i]); else printf("%d", tDep[i]);
        printf("])...\n");
        for(int i=3;i>0;i--){
			if (tDep[i] < tDep[i-1]){
				fprintf(stderr, "err: start(i) must be > start(i-1). Here %d < %d\n", tDep[i], tDep[i-1]);
				exit(1);
			}
		}
        printf("v=%d allOne (starting at #%llu)...\n", v, nbCfgDepart);
        singleOne(&jeu, v);

        // printf("recur...\n"); recur(0, v, &jeu);
        printf("Finished. nbCfg = %llu\n", nbCfg);
        clock_t stop = clock();
        printf("Finished in %6.3f s.\n", (double)(stop - start) / CLOCKS_PER_SEC);
        exit(0);
	}
	if (partialOn){
		int iv = jeuFromStr(&jeu, v, partial);
		printf("v=%d iv=%d ", v, iv); jeuDisp(&jeu); printf("\n");
		for(int i=1;i<=v;i++) tOcc[i]=0;
		for(int i=0;i<iv;i++) tOcc[jeu.cartes[i]]++;
		printf("tOcc=["); for(int i=1;i<=v;i++) printf("%d:%d ", i, tOcc[i]); printf("]\n");
		rFill(&jeu, v, iv, v*4);
		exit(1);
	}

    if (geneticOn){
		printf("loop genetic random configs (seed=%u nbCartes=%d roundLimit=%d)...\n", seed, nbCartes, roundLimit);
		for (int i = 0, k = 0; i < nbCartes; i += c) {
			for(int j = 0; j < c; j++) jeu.cartes[k++] = valeur;
			valeur++;
		}
		jeu.nb = nbCartes;
		jeuDisp(&jeu); printf("\n");
		printf("*** %d couleurs, %d valeurs ***\n", c, v);

		int populationSize = 32000;
		int nbParents = populationSize * parentsPercentage / 100;

		printf("populationSize=%d nbParents=%d mutationRatePer1000=%d/1000\n",  populationSize, nbParents, mutationRatePer1000);

		PEOPLE *population=(PEOPLE*)malloc(populationSize*sizeof(PEOPLE));
		// init population
		for(int i=0;i<populationSize;i++){
			jeuCopyTo(&jeu, &population[i].jeu);
			jeuShuffle(&population[i].jeu, 25);
			// printf("population[%d].jeu=",  i); jeuDisp(&population[i].jeu); printf("\n");
		}
		int nbCartesParJoueur = nbCartes / 2;
		time_t tDep = time(NULL), tNext = tDep + 2;
		// printf("tNext=%ld\n", tNext);
		int bestPrevious = 0, bestRepeat = 0, meteoritAfter=400, roundMaxLocal=0;
		for (uint generation=1;;generation++) {
			bool disp = false;
			if (roundMaxLocal == bestPrevious){
				bestRepeat ++;
				if (bestRepeat == meteoritAfter){
					// reset population
					fprintf(stderr, "BIG METEORIT after %d generations without progress...%c", bestRepeat, 13);
 					for(int i=0;i<populationSize;i++){
						jeuCopyTo(&jeu, &population[i].jeu);
						jeuShuffle(&population[i].jeu, nbCartesParJoueur);
					}
					bestRepeat = 0;
					if (meteoritAfter < 1000) meteoritAfter += 100;
				}
			}
			else{
				bestPrevious = roundMaxLocal;
				bestRepeat = 0;
			}

			time_t tCur = time(NULL);
			if (tCur >= tNext){ tNext = tCur + 2; disp = true; }
			// printf("tCur=%ld tNext=%ld delta=%ld disp=%d\n", tCur, tNext, tNext - tCur, disp?1:0);
			// if (disp) printf("gen#%u %ld %ld %c", generation, tCur, tNext, 13);
			// computes score for each person
			for(int i=0;i<populationSize;i++){
				distrib(&population[i].jeu);	// updates j1 and j2
				// printf("J1.nb=%d ", j1.nb); jeuDisp(&j1); printf(" "); printf("J2.nb=%d ", j2.nb); jeuDisp(&j2); printf("\n");
				TYP win = play(false);		// plays j1 and j2, updates iRound
				if (win == 0) {
					population[i].score = 0;	// printf("Deuce\n");
				}
				else if (win == 3) {
					population[i].score = 0;	// printf("Infinite loop\n");
				}
				else {
					// printf("*** winner : %d iRound=%d\n", win, iRound);
					population[i].score = iRound;
					// printf("iRound=%d J1=", iRound); jeuDisp(&j1); printf(" "); printf("J2.nb=%d ", j2.nb); jeuDisp(&j2); printf("\n");
					if (iRound > roundMaxLocal) {
						roundMaxLocal = iRound;
					}
					if (iRound > roundMax) {
						newBestScore(v, &population[i].jeu, iRound, generation);
					}
				}
			}
			// sort scores
			qsort(population, populationSize, sizeof(PEOPLE), cmpPeople);
			// display best
			if (disp) {
				fprintf(stderr, "#%u (noprogess=%d/%d) [", generation, bestRepeat, meteoritAfter);
				affDelai(tCur-tDep, stderr);
				fprintf(stderr, "] top10  =");
				for(int i=0;i<10;i++) fprintf(stderr, " %d", population[i].score);
				fprintf(stderr, "   %c", 13);
			}
			// keep 10 best
			int i = nbParents;
			// make new generations
			for(;i<populationSize;i++){
				// choose 2 parents, mix a part of p1 and the rest picked in p2 in order
				int iPere1 = i%nbParents, iPere2 = (iPere1+r1n(nbParents-1))%nbParents, nb1 = r1n(nbCartes-2) + 1;
				TYP *jcCur = population[i].jeu.cartes, *jc1 = population[iPere1].jeu.cartes;
				memcpy(jcCur, jc1, nb1 * sizeof(TYP));
				int tnb[16];
				for(int j=0;j<16;j++) tnb[j]=0;
				for(int j=0;j<nb1;j++) tnb[jc1[j]]++;
				TYP *jc2 = population[iPere2].jeu.cartes;
				for(int j=nb1, k=nb1;j<nbCartes;j++){
					while(tnb[jc2[k]]==4){ k++; k%=nbCartes; }	// deja 4 cartes, on passe a la suivante
					jcCur[j] = jc2[k];
					// printf(" - ajoute %d (nb=%d)\n", population[iPere2].jeu.cartes[k], tnb[population[iPere2].jeu.cartes[k]]);
					tnb[jc2[k]]++;
					k++; k%=nbCartes;
				}
/*
				// check for game integrity (4 cards)
				for(int j=0;j<nbCartes;j++) tnb[j]=0;
				for(int j=0;j<nbCartes;j++) tnb[jcCur[j]]++;
				for(int j=0;j<nbCartes;j++){
					int nb=tnb[jcCur[j]];
					if (nb>4){
						fprintf(stderr,"Oops, il y a %d cartes de valeur %d\n", nb, jcCur[j]);
						fprintf(stderr,"iPere1=%d iPere2=%d nb1=%d i=%d populationSize=%d\n", iPere1, iPere2, nb1, i, populationSize);
						printf("%4d", iPere1); jeuDisp(&population[iPere1].jeu); printf("\n");
						printf("%4d", iPere2); jeuDisp(&population[iPere2].jeu); printf("\n");
						printf("%4d", i); jeuDisp(&population[i].jeu); printf("\n");
						printf("iPere1     :"); for(int j=0;j<nb1;j++) printf(" %2d", jc1[j]); printf("\n");
						printf("Debut de i :"); for(int j=0;j<nb1;j++) printf(" %2d", jcCur[j]); printf("\n");
						for(int k=1;k<=nbCartes/4;k++) printf(" - tnb[%d]=%d\n", k, tnb[k]);
						exit(6);
					}
				}
*/
				/*
				// interchange x values of J1 with J2, ex : J1 : [4,1,2,1,4,5,3,2,5,5] J2 : [3,5,3,4,4,2,2,1,3,1]                                                          1,3                        2,5
				jeuCopyTo(&population[i%nbParents].jeu, &population[i].jeu);
				int iCut = r1n(nbCartesParJoueur-1) - 1;
				int iCutEnd = iCut + r1n(nbCartesParJoueur-iCut);
				for(int j=iCut;j<iCutEnd;j++){
					TYP tmp = population[i].jeu.cartes[j];
					population[i].jeu.cartes[j] = population[i].jeu.cartes[j+nbCartesParJoueur];
					population[i].jeu.cartes[j+nbCartesParJoueur] = tmp;
				}
				*/
			}
			// make mutations
			if (mutationRatePer1000){
				int nbMut = 0;
				for(i=1;i<populationSize;i++){	// no mutation on the best score
					// interchange 2 random values ex : J1 : [4,1,2,1,4,5,3,2,5,5] J2 : [3,5,3,4,4,2,2,1,3,1]
					//                                            5           2
					//  5% => 50%o
					/*
					int nbMutToDo = 0;
					// for(int j=0;j<nbCartes;j++) if ((RND%1000) <= mutationRatePer1000) nbMutToDo++;
					nbMutToDo = ((nbCartes * (RND%1000)) / 1000 ) * (mutationRatePer1000 / 1000);
					int r=RND%101;
					nbMutToDo = (mutationRatePer1000 * nbCartes * r / 10 / 100 / 100); printf("r=%d nbMutToDo=%d\n", r, nbMutToDo);
					exit(9);
					for(int j=0;j<nbMutToDo;j++){
						nbMut++;
						int i1 = RND % nbCartes, i2 = (i1 + RND) % nbCartes;
						swapCartes(&population[i].jeu, i1, i2);
					}
					*/
					for(int j=0;j<nbCartes;j++){
						int r = RND;
						if ( (r%1000) < mutationRatePer1000){
							nbMut++;
							// swapCartes(&population[i].jeu, j, (j + 1 + (r % (nbCartes - 1))) % nbCartes);
							swapCartes(&population[i].jeu, j, r % nbCartes);
						}
					}
				}
				/*
				fprintf(stderr, "nbMut=%d/%d=%.1lf/indiv%% %d/%d=%.1lf/gene%% "
					, nbMut, populationSize, 100.0*(double)nbMut/(double)populationSize
					, nbMut, populationSize*nbCartes, 100.0*(double)nbMut/(double)(populationSize*nbCartes)
					);
				*/
			}
			// 2-10 : 1 mutation from best
			for(i=1;i<5*populationSize/100;i++){
				jeuCopyTo(&population[0].jeu, &population[i].jeu);
				// interchange 2 random values ex : J1 : [4,1,2,1,4,5,3,2,5,5] J2 : [3,5,3,4,4,2,2,1,3,1]
				//                                          2           1
				for(;;){
					int i1 = RND % nbCartes, i2 = (i1 + RND) % nbCartes;
					TYP tmp = population[i].jeu.cartes[i1];
					if (population[i].jeu.cartes[i2] == tmp) continue;
					swapCartes(&population[i].jeu, i1, i2);
					break;
				}
			}
			// small part totally random
			for(i=95*populationSize/100;i<populationSize;i++){
				jeuShuffle(&population[i].jeu, 25);
			}

			// printf("#%u nbMut=%d/%d   %c", generation, nbMut, populationSize, 13);
			// break;
		}
		exit(0);
	}

    if (!randomOn) {
        // loop all config
        printf("loop all configs (roundLimit=%d)...\n", roundLimit);
        start = clock();
        recurInit(v);
        if (v > 4) dispInfinite = true;
        // jeuDisp(&jeu); printf("\n");

        // printf("allOne...\n"); allOne(&jeu, v);

        printf("v=%d allOne (starting at #%llu)...\n", v, nbCfgDepart);
        allOne(&jeu, v);

        // printf("recur...\n"); recur(0, v, &jeu);
        printf("Finished. nbCfg = %llu\n", nbCfg);
        clock_t stop = clock();
        printf("Finished in %6.3f s.\n", (double)(stop - start) / CLOCKS_PER_SEC);
        exit(0);
    }

    if (v >= 5) dispInfinite = true;

	printf("loop random configs (method=%s seed=%u roundLimit=%d)...\n", randomMethod, seed, roundLimit);
    for (int i = 0, k = 0; i < nbCartes; i += c) {
        for(int j = 0; j < c; j++) jeu.cartes[k++] = valeur;
        valeur++;
    }
    jeu.nb = nbCartes;
    jeuDisp(&jeu); printf("\n");
    printf("*** %d couleurs, %d valeurs. Previous record : %d ***\n", c, v, record);
    bool iniCfg = strcmp(cfg, "nope");

    if (!strcmp(randomMethod, "mc2")){
		// essaie x solutions au hasard avec c1=1,2,3,4 garde la meilleure, et continue avec c2, c3...
		JEU jeuSave, jeuBest;
		printf("jeu="); jeuDisp(&jeu); printf("\n");
		jeuCopyTo(&jeu, &jeuSave);
		int resultBest = game(jeu, 1, v);
		static const int KMAX = 100000;
		time_t tDepart = time(NULL), tNow;
		char speedStr[256];
		strcpy(speedStr, "---");
		double attemptTot = 0.0;
		int iDep = 0;
		if (strlen(partial) != 0){
			jDep.nb = v * 4;
			int iv = jeuFromStr(&jDep, v, partial);
			printf("v=%d iv=%d jDep=", v, iv); jeuDisp(&jDep); printf("\n");
			for(int i=0;i<iv;i++){
				if (jeu.cartes[i] == jDep.cartes[i]) continue;
				int iSwap = -1;
				for(int k=i+1;k<nbCartes;k++) if (jeu.cartes[k]==jDep.cartes[i]){ iSwap = k; swapCartes(&jeu, i, k); break; }
				if (iSwap < 0){
					fprintf(stderr, "Error : can't find free card '%d' to put in place %d\n", jDep.cartes[i], i);
					exit(3);
				}
			}
			iDep = iv;
			jeuCopyTo(&jeu, &jDep);
			printf("v=%d iv=%d jDep=", v, iv); jeuDisp(&jDep); printf("\n");
		}
		for (uint attempt=1;;attempt++) {
			if (iDep == 0) jeuShuffle(&jeu, 24); else jeuCopyTo(&jDep, &jeu);
			// printf("pass #%8u jeu=", attempt); jeuDisp(&jeu); printf("  %c", 13);
			fprintf(stderr, "#%u ", attempt);
			jeuDispFile(stderr, &jeu);
			tNow = time(NULL);
			if (tNow > tDepart){
				dispDelai(stderr, tDepart, tNow);
				sprintf(speedStr, " %.2lf M/s", attemptTot/(double)(tNow-tDepart)/1000000.0);
				fprintf(stderr, "%s", speedStr);
			}
			fprintf(stderr, " %c", 13);

			for(int i=iDep;i<nbCartes-1;i++){
				// fprintf(stderr, "pass %u %d/%d  %c", attempt, i, nbCartes-1, 13);
				int jBest = 0;
				for(int j=1;j<=v;j++){
					int iSwap = -1;
					if (jeu.cartes[i] != j){
						for(int k=i+1;k<nbCartes;k++) if (jeu.cartes[k]==j){ iSwap = k; swapCartes(&jeu, i, k); break; }
						if (iSwap < 0) continue;
					}
					jeuCopyTo(&jeu, &jeuSave);
					for(int k=0;k<KMAX;k++){	// run 10000 random games
						// shuffle cards from i+1 to nbCartes
						for(int r=0;r<8;r++) swapCartes(&jeu, i+1+RND%(nbCartes-i-1), i+1+RND%(nbCartes-i-1));
						int result = game(jeu, attempt, v);
						if (result > resultBest){
							resultBest = result;
							jBest = j;
							jeuCopyTo(&jeu, &jeuBest);
							newBestScore(v, &jeu, result, attempt, speedStr);	// display best and eventually save record to disk
						}
					}
					attemptTot += (double)KMAX;
					jeuCopyTo(&jeuSave, &jeu);	// restore
					if (iSwap >= 0) swapCartes(&jeu, i, iSwap);	// unswap
				}
				if (jBest==0){
					// printf("i=%d/%d : rien de mieux...\n", i, nbCartes-1);
				}
				else{
					// printf("jeu="); jeuDisp(&jeu); printf("\n");
					// printf("play jBest=%d at pos %d\n", jBest, i);
					/*
					if (jeu.cartes[i] != jBest){
						for(int k=i+1;k<nbCartes;k++) if (jeu.cartes[k]==jBest){ swapCartes(&jeu, i, k); break; }
					}
					*/
					jeuCopyTo(&jeuBest, &jeu);
					printf("jeu="); jeuDisp(&jeu, i); printf("\n");
				}
			}
		}
		exit(0);
	}

    if (!strcmp(randomMethod, "montecarlo")){
		// essaie des solutions au hasard, et pour chacune essaie toutes les combinaisons de swap, et recommence avec la meilleure, sinon remelange tout
		jeuShuffle(&jeu, 24);
		JEU jeuSave, jeuBest, jeuBestEver;
		time_t tDepart = time(NULL);
		for (uint attempt=1;;attempt++) {
			if ((attempt%10000)==0){
				fprintf(stderr, "#%u ", attempt);
				dispDelai(stderr, tDepart, time(NULL));
				fprintf(stderr, " %c",13);
			}
			jeuCopyTo(&jeu, &jeuSave);
			int resultLocal = game(jeu, attempt, v);
			bool betterFound = false;
			for(int i=0;i<nbCartes-1;i++){
			   for(int j=i+1;j<nbCartes;j++){
				   if (jeu.cartes[i] == jeu.cartes[j]) continue;
					TYP tmp = jeu.cartes[i];	jeu.cartes[i] = jeu.cartes[j]; jeu.cartes[j] = tmp;	// swap 2 cards
					int result = game(jeu, attempt, v);
					if (result > resultLocal){
						resultLocal = result;
						jeuCopyTo(&jeu, &jeuBest);
						betterFound = true;
					}
					tmp = jeu.cartes[i]; jeu.cartes[i] = jeu.cartes[j]; jeu.cartes[j] = tmp;		// restore
				}
			}
			if (betterFound){
				/*
				if (attempt<10){
					printf("better found : score=%d ", resultLocal);
					distrib(&jeuBest);
					printf(" J1 : ");	jeuDisp(&j1);
					printf(" J2 : ");	jeuDisp(&j2);
					printf("\n");
                }
                */
				jeuCopyTo(&jeuBest, &jeu);
			}
			else{
				// if (attempt<10){ printf("shuffle.\n"); }
				jeuShuffle(&jeu, 24);
			}
		}
		exit(0);
	}

    time_t t0=time(NULL);
    for (uint attempt=1;;attempt++) {
        if ((attempt%1000000)==0) fprintf(stderr, "#%u %c", attempt, 13);
        jeuShuffle(&jeu, 16);
        if (iniCfg) {
            int nbCfg = cleanupCfg(cfg, &jeu);
            if (nbCfg != nbCartes) {
                fprintf(stderr, "bad cfg length : %d insteads of %d\n", nbCfg, nbCartes);
                return 1;
            }
            // jeu.nb = nbCfg;
            // for (int i = 0; i < jeu.nb; i++) jeu.cartes[i] = cfg[i] - '0';
            jeuDisp(&jeu); printf("\n");
        }
        // jeuCopyTo(&jeu, &jeuCopy);
        // jeuDisp(&jeu); printf("\n");
        distrib(&jeu);					// updates j1 and j2
        TYP win = play(false);		// plays j1 and j2, updates iRound
        if (win == 0) {
            // printf("Deuce\n");
        }
        else if (win == 3) {
			if (iRound < loopSizeMin){
				loopSizeMin = iRound;
				// printf("******************************\n");
				printf("*** %d : new min loop size %d", jeu.nb/4, loopSizeMin);
				// printf(" for jeu="); jeuDisp(&jeu); printf("\n");
				printf(" [%d", jeu.cartes[0]); for(int i=1;i<nbCartes/2;i++) printf(",%d", jeu.cartes[i]);
				printf("] [%d", jeu.cartes[nbCartes/2]); for(int i=nbCartes/2+1;i<nbCartes;i++) printf(",%d", jeu.cartes[i]);
				printf("]\n");
				if (loopSizeMin < loopMinRecord){
					char buf[2048], word[512];
					sprintf(buf, "[%d", jeu.cartes[0]); for(int i=1;i<nbCartes/2;i++){ sprintf(word, ",%d", jeu.cartes[i]); strcat(buf, word); }
					sprintf(word, "] [%d", jeu.cartes[nbCartes/2]); strcat(buf, word); for(int i=nbCartes/2+1;i<nbCartes;i++){ sprintf(word, ",%d", jeu.cartes[i]); strcat(buf, word); }
					strcat(buf, "]");
					// strcat(buf, "\n");
					recordSave(jeu.nb/4, loopSizeMin, "loop", buf);
				}
				// printf("******************************\n");
				// printf("loop size %d detected\n", iRound);
				// for(int i=0;i<iLoopMax;i++){ printf("%4d", i); jeuDisp(&loop[i].j1); printf(" "); jeuDisp(&loop[i].j2); printf(" "); jeuDisp(&loop[i].c); printf("\n"); }
				// int k = iRound-1; printf("%4d", k); jeuDisp(&loop[k].j1); printf(" "); jeuDisp(&loop[k].j2); printf(" "); jeuDisp(&loop[k].c); printf("\n");
				// k -= iRound; printf("%4d", k); jeuDisp(&loop[k].j1); printf(" "); jeuDisp(&loop[k].j2); printf(" "); jeuDisp(&loop[k].c); printf("\n");
				// exit(1);
			}
            // printf("Infinite loop iRound=%d", iRound);
        }
        else {
            // printf("*** winner : %d iRound=%d\n", win, iRound);
            if (iRound > roundMax) {
                roundMax = iRound;
                // jeuDisp(&jeu); printf(" :: ");
                distrib(&jeu);
                printf("J1 : ");        jeuDisp(&j1);
                printf(" J2 : ");       jeuDisp(&j2);
                printf(" : new nbRoundMax=%d/%d", roundMax, roundLimit);
                printf(" after %u attempts", attempt);
                time_t t1 = time(NULL);
                printf(" (%ld')", (t1-t0));
                if (t1==t0) printf("---M/s"); else printf(" (%.2lfM/s)", (double)attempt/(double)(t1-t0));
                printf("\n");
                if (roundMax > record) {
                    roundLimit = roundMax + 100;	// update roundLimit
                    record = roundMax;
                    char buf[512], word[256];
                    jeuDispStr(&j1, buf);		strcat(buf, " ");
                    jeuDispStr(&j2, word);		strcat(buf, word);
                    sprintf(word, " after %u attempts", attempt);	strcat(buf, word);
                    recordSave((uint)v, roundMax, "len", buf);
                }
            }
        }
        if (iniCfg) {
            break;
        }
    }
    printf("Finished.\n");
}

