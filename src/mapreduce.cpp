/*
 * points.h
 *
 *  Created on: Jul 16, 2013
 *      Author: ricardo
 */

#define MESSAGE_BUFFER_SIZE 100000

#define USE_ALL_FILES_LOCALLY true
#define MAX_DIMENSIONS 5

#include <vector>
#include <math.h>
#include <string>

#include <iostream>
#include "mpi.h" 
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>

#include <dirent.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include <new>


namespace uff{
    class Point;
    class PointPair;
    class Serialization;
    class Calc;
    struct KeyPair;
}

//// Point ------------------------------

class uff::Point {
    private:
        int dimension;
        double *values;
        long index;
        long pid;

    public:
        Point(int);
        Point();
        void setValue(int, double);
        double getValue(int);
        void setIndex(long i);
        long getIndex();
        void setPid(long p);
        long getPid();
        std::string serialize();
        int getDimension();
        bool equals(Point* p);
};

class uff::PointPair {
    public:
        uff::Point *pointA;
        uff::Point *pointB;
        double similarity;

        PointPair(uff::Point* pA, uff::Point* pB, double sim);
        ~PointPair();

        std::string serialize();
};

struct uff::KeyPair{
    int value1;
    int value2;

    std::string serialize();
};


class uff::Serialization{
    public:
        static uff::Point* unserializePoint(std::string strPoint);
        static uff::PointPair* unserializePointPair(std::string strPointPair);
        static uff::KeyPair unserializeKey(std::string strKey);
};


class uff::Calc{
    public:
        static double calculateSimilarity(uff::Point* pA, uff::Point* pB);
};


/*
 * points.cpp
 *
 *  Created on: Jul 17, 2013
 *      Author: ricardo
 */

uff::Point::Point(int dim){
    dimension = dim;
    values = (double *) malloc(dimension * sizeof(double));
    pid = -1;
    index = -1;
    for (int i=0; i<dimension; i++){
        values[i] = 0.0;
    }
}
uff::Point::Point(){
    pid = -1;
    index = -1;
    dimension = 0;
    values = NULL;
}

void uff::Point::setValue(int i, double value){
    values[i] = value;
}

double uff::Point::getValue(int i){
    return values[i];
}

void uff::Point::setIndex(long i){
    index = i;
}
long uff::Point::getIndex(){
    return index;
}
void uff::Point::setPid(long p){
    pid = p;
}
long uff::Point::getPid(){
    return pid;
}

std::string uff::Point::serialize(){
    std::string sb = "";
    char strValue[256];
    if (index >= 0){
        sprintf(strValue, "%li", index);
        sb.append(strValue);
        if (pid >= 0){
            sb.append("/");
            sprintf(strValue, "%li", pid);
            sb.append(strValue);
        }
        sb.append(":");
    }
    sb.append("(");
    std::string comma = "";
    for (int i=0; i<dimension; i++){
        sprintf(strValue, "%g", values[i]);
        sb.append(comma);
        sb.append(strValue);
        comma = ", ";
    }
    sb.append(")");
    return sb;
}

int uff::Point::getDimension(){
    return dimension;
}

bool uff::Point::equals(uff::Point* point){
    bool result = true;
    result = result && (this->dimension == point->dimension);
    for(int i=0; i<dimension; i++){
        result = result && (this->values[i] == point->getValue(i));
    }
    return result;
}


/**
 * Format: index/pid:(x1, y1, ...)
 */
uff::Point* uff::Serialization::unserializePoint(std::string strPoint){
    int charPos = strPoint.find(":");
    int index = 0;
    int pid = 0;
    if (charPos != std::string::npos){
        std::string strIndex = strPoint.substr(0, charPos);
        sscanf(strIndex.c_str(), "%i/%i", &index, &pid);
    }

    charPos = strPoint.find("(");
    if (charPos == std::string::npos){
        return NULL;
    }
    strPoint = strPoint.substr(charPos+1, strPoint.length() - charPos);
    charPos = strPoint.find(")");
    if (charPos == std::string::npos){
        return NULL;
    }
    strPoint = strPoint.substr(0, charPos);
    int dim = 1;
    for(unsigned int i=0; i<strPoint.size(); i++){
        if (strPoint[i]==','){
            dim++;
        }
    }
    uff::Point* pt = new uff::Point(dim);
    pt->setIndex(index);
    pt->setPid(pid);

    double number;
    std::string strNumber;
    int pos;
    for (int i=0; i<dim; i++){
        pos=strPoint.find(", ");
        if (pos==std::string::npos){
            pos=strPoint.length();
        }
        strNumber = strPoint.substr(0, pos);
        number = atof(strNumber.c_str());
        strPoint = strPoint.erase(0, pos + 1);
        pt->setValue(i, number);
    }
    return pt;
}

// End Point ---------------------------




//// PointPair ---------------------------

uff::PointPair::PointPair(uff::Point* pA, uff::Point* pB, double sim){
    pointA = pA;
    pointB = pB;
    similarity = sim;
}

uff::PointPair::~PointPair(){

}

std::string uff::PointPair::serialize(){
    std::string sb = "";
    sb.append("{");
    sb.append(pointA->serialize());
    sb.append(";");
    sb.append(pointB->serialize());
    sb.append(" s: ");
    char strValue[256];
    sprintf(strValue, "%g", similarity);
    sb.append(strValue);
    sb.append("}");
    return sb;
}

/*
 * Format: {index1/pid1:(p1, p2);index2/pid2:(p1, p2) s: similarity}
 */
uff::PointPair* uff::Serialization::unserializePointPair(std::string strPointPair){
    int charPos = strPointPair.find("s:");
    if (charPos == std::string::npos){
        return NULL;
    }
    std::string strSimilarity = strPointPair.substr(charPos + 3, strPointPair.length() - charPos - 4);
    double similarity = atof(strSimilarity.c_str());
//    std::cout << "M" << strPointPair << " - *" << strSimilarity << "*" <<  similarity << std::endl;
    sscanf(strPointPair.c_str(), "s: %lg}", &similarity);
    int startPos = strPointPair.find("{");
    if (startPos == std::string::npos){
        startPos = 0;
    }else{
        startPos = 1;
    }
    std::string strPoints = strPointPair.substr(startPos, charPos - startPos);
    charPos = strPoints.find(";");
    std::string strPoint1 =  strPoints.substr(0, charPos);
    std::string strPoint2 =  strPoints.substr(charPos + 1, strPoints.length() - charPos);

    uff::Point* point1 = uff::Serialization::unserializePoint(strPoint1);
    uff::Point* point2 = uff::Serialization::unserializePoint(strPoint2);

    char strValue[256];
    sprintf(strValue, "%g", similarity);

//    std::cout << "X:" << strPointPair << " " << strPoints << " " << strValue << std::endl;
    uff::PointPair* pair = new uff::PointPair(point1, point2, similarity);
    return pair;
}

// End PointPair -----------------------

std::string uff::KeyPair::serialize(){
    char str[500];
    sprintf(str, "(%i, %i)", value1, value2);
    return std::string(str);
}


// Domain functions ------------------------------------

uff::KeyPair uff::Serialization::unserializeKey(std::string strKey){
    uff::KeyPair pair;
    sscanf(strKey.c_str(),"[%d,%d]", &pair.value1, &pair.value2);
    return pair;
}

double uff::Calc::calculateSimilarity(uff::Point* pointA, uff::Point* pointB){
    if (pointA->getDimension() != pointB->getDimension()){
        return -1;
    }
    double squareSum = 0;
    for(int i=0; i<pointA->getDimension(); i++){
        squareSum += pow(pointA->getValue(i) - pointB->getValue(i), 2);
    }
    double result = sqrt( squareSum );
    return result;
}


#define START_MAPPER_MESSAGE 1
#define MAP_MESSAGE 10
#define RUN_REDUCE_MESSAGE 20
#define REDUCE_MESSAGE 30
#define HOSTNAME_SIZE 255

namespace uff {
class TopKSimilarityMapRedApp;
}

class uff::TopKSimilarityMapRedApp {
private:
    int k, m;
    int numReducers, numMappers;
    char hostname[HOSTNAME_SIZE];

    int localRank, localSize;
    int mapperStartRank, reducerStartRank;
    int mpiRank, mpiSize;
    std::vector<std::string> filesToMap;
    std::string* messages;
    std::vector<uff::Point*> pointList;
    std::vector<uff::PointPair*>* sortedPairs;

    uff::KeyPair* reduceKey;

    int bufferSize;

    std::string localInputFolder;
    std::string outpuFolder;

public:
    TopKSimilarityMapRedApp(int size, int rank);
    void run(int k, int m, int numFiles, std::string folder, std::string outputFolder);
    void startMappers();
    void waitForMaster();
    void discoverLocalRank();
    void selectInputFiles();
    void runMappers();
    void registerMessageToReduce(int reduceIndex, std::string key,
            std::string value);
    void sendAllMessagesToReducers();
    void sendOneMessageToReducer(int reduceIndex, std::string message);
    void receiveMessagesFromMappers();
    void extractKeyValuesFromMessage(char* message);
    void extractPointPairsFromMessage(char* message, std::vector<uff::PointPair*>* receivedSortedPairs);
    void mergeSortedArrays(std::vector<uff::PointPair*>* receivedSortedPairs);
    void runReducers();
    void runFinalBinaryReduce();
    void writeResultToFile();
    void receiveFinalResultFromReduce();
    void sendResultToMaster();

    void checkForSendErrors(int ierr, std::string msg);
    void log(std::string msg);

    std::string* serializePointList(std::vector<uff::Point*>* pointList);
    std::string* serializePairList(std::vector<uff::PointPair*>* pairList);

    void pointToPairMapper(long index, Point* p);
    void topKDistanceReducer(uff::KeyPair* key,
            std::vector<uff::Point*> *points);
};

uff::TopKSimilarityMapRedApp::TopKSimilarityMapRedApp(int size, int rank) {
    this->mpiSize = size;
    this->mpiRank = rank;
    this->k = 0;
    this->m = 0;
    this->numMappers = 0;
    this->numReducers = 0;
    this->reduceKey = NULL;
    this->messages = NULL;
    this->localRank = 0;
    this->localSize = 0;
    this->sortedPairs = NULL;
    this->mapperStartRank = 1;
    this->reducerStartRank = 1;
    this->bufferSize = 0;
    gethostname(this->hostname, HOSTNAME_SIZE);
}

void uff::TopKSimilarityMapRedApp::run(int k, int m, int numFiles, std::string inputFolder, std::string output) {
    this->k = k;
    this->m = m;

    this->numReducers = (m * (m + 1)) / 2;
    this->numMappers = numFiles;

    int allProcess=1 + (this->numReducers + this->numMappers);
    if ( allProcess > mpiSize){
        if (mpiRank == 0){
            std::cerr << "Number of process is wrong. It should be " << allProcess << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 2);
        return;
    }


    this->mapperStartRank = 1; // Master process
    this->reducerStartRank = 1 + this->numMappers;
    this->messages = new std::string[numReducers];
    this->localInputFolder = inputFolder;
    this->outpuFolder = output;

    this->bufferSize = 25 *  10 * k;

    // Master ----------------------------------------------------------------
    if (mpiRank == 0){
        this->log("Root process starting ");
        this->startMappers();
try{
        this->receiveFinalResultFromReduce();
} catch (std::bad_alloc& ba) {
    std::cerr << mpiRank << " master receiveFinalResultFromReduce bad_alloc caught: " << ba.what() << std::endl;
    return;
}
        std::cout << mpiRank << " 8.5 - Write to file." << std::endl;
        this->writeResultToFile();

    }
    // Mapper code -----------------------------------------------------------
    else if (mpiRank < reducerStartRank) {

        this->log("0.0 - Waiting for master to start...");
        this->waitForMaster();

        // Discover local rank (unique rank per machine)
        this->log("1.0 - Discover local rank...");
        this->discoverLocalRank();

        // Iterate over files and select what files the current map must use as input
        this->log("2.0 - Selecting input files...");
        this->selectInputFiles();

        // Run map for each line in the selected files.
        this->log("3.0 - Run mappers ...");
        this->runMappers();

        // Send and receive mapper messages
        this->log("4.0 - Send all map messages...");
        this->sendAllMessagesToReducers();

        this->log("4.5 - MAP FINALIZING...");

    }
    // Reduce code ------------------------------------------------------------
    else if (mpiRank >= reducerStartRank
            && mpiRank < reducerStartRank + numReducers){

try{
        // Receive messages from mappers
        this->log("5.0 - Waiting for messages from mappers...");
        this->receiveMessagesFromMappers();
} catch (std::bad_alloc& ba) {
    std::cerr << mpiRank << " receiveMessagesFromMappers bad_alloc caught: " << ba.what() << std::endl;
    return;
}
try{
        // Run reducer for all received points
        this->log("6.0 - Run reducers");
        this->runReducers();
} catch (std::bad_alloc& ba) {
    std::cerr << mpiRank << " runReducers bad_alloc caught: " << ba.what() << std::endl;
    return;
}
try{
        // Run binary reduce to merge all top k pairs
        this->log("7.0 - Run final binary reduce...");
        this->runFinalBinaryReduce();
} catch (std::bad_alloc& ba) {
    std::cerr << mpiRank << " runFinalBinaryReduce bad_alloc caught: " << ba.what() << std::endl;
    return;
}
try{
        if (mpiRank == reducerStartRank){
            this->sendResultToMaster();
        }
} catch (std::bad_alloc& ba) {
    std::cerr << mpiRank << " sendResultToMaster bad_alloc caught: " << ba.what() << std::endl;
    return;
}
    }
}

void uff::TopKSimilarityMapRedApp::startMappers() {
    int intMessage=1;
    for (int i = 0; i < numMappers; i++) {
        int ierr=MPI_Send(&intMessage, 1, MPI_INT, mapperStartRank + i, START_MAPPER_MESSAGE, MPI_COMM_WORLD);
        this->checkForSendErrors(ierr, "uff::TopKSimilarityMapRedApp::startMappers");
    }
}

void uff::TopKSimilarityMapRedApp::waitForMaster() {
    MPI_Status status;
    int i;
    MPI_Recv(&i, 1, MPI_INT, 0, START_MAPPER_MESSAGE, MPI_COMM_WORLD, &status);
}

void uff::TopKSimilarityMapRedApp::discoverLocalRank() {
    char remoteHostname[HOSTNAME_SIZE];
    MPI_Status status;
    localSize = 1; // current process counts as 1
    localRank = 0;
    for (int i = mapperStartRank; i < (mapperStartRank + numMappers); i++) {
        if (i == mpiRank)
            continue;
        if (i < mpiRank) {
            MPI_Send(this->hostname, HOSTNAME_SIZE, MPI_CHAR, i, 1, MPI_COMM_WORLD);
            MPI_Recv(remoteHostname, HOSTNAME_SIZE, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(remoteHostname, HOSTNAME_SIZE, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
            MPI_Send(this->hostname, HOSTNAME_SIZE, MPI_CHAR, i, 1, MPI_COMM_WORLD);
        }
        if (strcmp(remoteHostname, this->hostname) == 0) {
            localSize++;
            if (status.MPI_SOURCE < this->mpiRank) {
                localRank++;
            }
        }
    }
    std::cout << mpiRank << " 1.1 " << hostname << " localRank: " << localRank << " size: " << localSize << std::endl;
}

void uff::TopKSimilarityMapRedApp::selectInputFiles() {
    DIR *dir = opendir(localInputFolder.c_str());
    if (dir == NULL) {
        std::cerr << "Can't open folder '" << localInputFolder << "'"
                    << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }
    struct dirent *ent;
    unsigned long fileIndex = 0;
    std::vector<std::string>::iterator it;

    std::vector<std::string> sortedFiles;
    while ((ent = readdir(dir)) != NULL) {
        std::string fullFileName = localInputFolder + "/" + ent->d_name;
        DIR *dir2 = opendir(fullFileName.c_str());
        if (dir2 != NULL || ent->d_name[0] == '.') {
            continue;
        }
        // Adicionar a lista de arquivos ordenada
        if (sortedFiles.size() == 0) {
            sortedFiles.push_back(ent->d_name);
        } else {
            bool inserted = false;
            for (it = sortedFiles.begin(); it < sortedFiles.end(); it++) {
                if (ent->d_name < (*it)) {
                    sortedFiles.insert(it, ent->d_name);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) sortedFiles.push_back(ent->d_name);
        }
    }
    for (fileIndex=0; fileIndex<sortedFiles.size(); fileIndex++) {
        // Does the node have all files?
        if (USE_ALL_FILES_LOCALLY){
            if (fileIndex % numMappers == (mpiRank - mapperStartRank)) {
                std::cout << "FILE: " << mpiRank << " " << sortedFiles[fileIndex] << std::endl;
                filesToMap.push_back(sortedFiles[fileIndex]);
            }
        } else {
            // Distribute files between local mappers
           if (fileIndex % localSize == localRank) {
               filesToMap.push_back(sortedFiles[fileIndex]);
           }
        }
    }
    closedir(dir);
}

void uff::TopKSimilarityMapRedApp::runMappers() {
    std::string line;
    std::string strIndex;
    Point* pt;

    // For each file for the mapper
    for (std::vector<std::string>::iterator it = filesToMap.begin();
            it != filesToMap.end(); ++it) {
        std::string fileName = *it;
        std::string fullFileName = localInputFolder + "/" + fileName;

        std::cout << mpiRank << " 3.1 " << fullFileName << std::endl;

        std::ifstream pointsFile(fullFileName.c_str());
        // Read file line by line
        if (pointsFile.is_open()) {
            while (pointsFile.good()) {
                getline(pointsFile, line);

                int charPos = line.find(":");
                if (charPos == std::string::npos) {
                    continue;
                }
                pt = uff::Serialization::unserializePoint(line);
                if (pt->getDimension() > 0) {
                    pointToPairMapper(pt->getIndex(), pt);
                }
            }
            pointsFile.close();
        }
    }
}

void uff::TopKSimilarityMapRedApp::registerMessageToReduce(int reduceIndex,
        std::string key, std::string value) {
//    cout << mpiRank << " reduce: " << reduceIndex << " " << key << " " << value << endl;
    if (messages[reduceIndex].length() == 0) {
        messages[reduceIndex].append(key);
    } else {
        messages[reduceIndex].append("|");
    }
    messages[reduceIndex].append(value);
    if (messages[reduceIndex].length() >= (MESSAGE_BUFFER_SIZE - value.length())){
        sendOneMessageToReducer(reduceIndex, messages[reduceIndex]);
        messages[reduceIndex]="";
    }
}


void uff::TopKSimilarityMapRedApp::sendAllMessagesToReducers() {
    char runReduceMessage[] = "RUN_REDUCE_MESSAGE";
    // Send one message, with all points, for each reduce
    for (int i = 0; i < numReducers; i++) {
        if (messages[i].length() > MESSAGE_BUFFER_SIZE) {
            std::cerr << "Buffer size exceed." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
            exit(2);
        }
        if (messages[i].length() > 0) {
            sendOneMessageToReducer(i, messages[i]);
        }
        MPI_Send(runReduceMessage, strlen(runReduceMessage) + 1, MPI_CHAR, (i + reducerStartRank), RUN_REDUCE_MESSAGE, MPI_COMM_WORLD);
    }
    std::cout << mpiRank << " 4.9 message sent" << std::endl;
}

void uff::TopKSimilarityMapRedApp::sendOneMessageToReducer(int reduceIndex, std::string message) {
    std::cout << mpiRank << " 4.1 sending to " << reduceIndex << " size " << message.length() << std::endl;
    char messageBuffer[MESSAGE_BUFFER_SIZE];
    strcpy(messageBuffer, message.c_str());
    int ierr = MPI_Send(messageBuffer, MESSAGE_BUFFER_SIZE, MPI_CHAR, reduceIndex + reducerStartRank, MAP_MESSAGE, MPI_COMM_WORLD);
    this->checkForSendErrors(ierr, "uff::TopKSimilarityMapRedApp::sendOneMessageToReducer");
}

void uff::TopKSimilarityMapRedApp::receiveMessagesFromMappers(){
    char *messageBuffer = new char[MESSAGE_BUFFER_SIZE];
    MPI_Status status;
    // Receive messages from all mappers
    int runMessagesCount = 0;
    while (runMessagesCount < numMappers){
        MPI_Recv(messageBuffer, MESSAGE_BUFFER_SIZE, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == MAP_MESSAGE){
            std::cout  << mpiRank << " 5.1 received MAP_MESSAGE from " << status.MPI_SOURCE << std::endl;
            this->extractKeyValuesFromMessage(messageBuffer);
        } else {
            std::cout  << mpiRank << " 5.1 received RUN_REDUCE_MESSAGE from " << status.MPI_SOURCE << std::endl;
            runMessagesCount++;
        }
    }
}

void uff::TopKSimilarityMapRedApp::extractKeyValuesFromMessage(char* message) {

    KeyPair key;
    int parseStatus = sscanf(message, "[%i, %i", &key.value1, &key.value2);
    if (parseStatus <= 0) {
        return;
    }
    // First message?
    if (this->reduceKey == NULL) {
        this->reduceKey = new KeyPair();
        this->reduceKey->value1 = key.value1;
        this->reduceKey->value2 = key.value2;

        // Current Key must match previous one
    } else if (this->reduceKey->value1 != key.value1
            || this->reduceKey->value2 != key.value2) {
        std::cout << "Reduce keys don't match " << this->reduceKey->serialize()
                << " new: " << key.serialize() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
        exit(2);
    }

    std::string strPoints = std::string(message);
    int charPos = strPoints.find("]");
    strPoints = strPoints.substr(charPos + 1, strPoints.length() - charPos);

    std::string strPoint;
    uff::Point* point;
    std::vector<Point*> list;
    while (strPoints.length() > 0) {
        charPos = strPoints.find('|');
        if (charPos == std::string::npos) {
            charPos = strPoints.length();
        }
        strPoint = strPoints.substr(0, charPos);
        point = uff::Serialization::unserializePoint(strPoint);

        this->pointList.push_back(point);

        strPoints = strPoints.erase(0, charPos + 1);
    }
}

void uff::TopKSimilarityMapRedApp::extractPointPairsFromMessage(char* message, std::vector<uff::PointPair*>* receivedSortedPairs) {
    std::string strMessage = std::string(message);
    receivedSortedPairs->clear();
    int charPos;
    std::string strPair;
//    std::cout << mpiRank << " - " << message << std::endl;
    PointPair* pair;
    while(strMessage.length() > 0){
        charPos = strMessage.find("|");
        if (charPos == std::string::npos) {
            charPos = strMessage.length();
        }
        strPair = strMessage.substr(0, charPos);
        pair = Serialization::unserializePointPair(strPair);
        receivedSortedPairs->push_back(pair);
//        std::cout << mpiRank << " - " << strPair << std::endl;
//        std::cout << mpiRank << " Pair: " << pair->serialize() << std::endl;
        strMessage = strMessage.erase(0, charPos + 1);
    }
}


void uff::TopKSimilarityMapRedApp::runReducers() {
    if (this->reduceKey != NULL && this->reduceKey->value1 > 0
            && this->reduceKey->value2 > 0) {
        this->topKDistanceReducer(reduceKey, &pointList);
    }
//    std::string* message = serializePairList(sortedPairs);;
//    std::cout << "** " << mpiRank << " " << message->c_str() << std::endl;
}

void uff::TopKSimilarityMapRedApp::runFinalBinaryReduce() {
    int numRounds = ceil(log2(this->numReducers));
    int previousPivot = 1;
    int indexPivot = 2;
    int sender, receiver;
    MPI_Status status;
    int reducerRank = mpiRank - reducerStartRank;

    char* messageBuffer;

    std::vector<uff::PointPair*>* receivedSortedPairs = new std::vector<uff::PointPair*>();
    std::string* message;


    for (int i = 1; i <= numRounds; i++) {
        if (mpiRank == reducerStartRank) std::cout << mpiRank << " 7.1 - Round: " << i << std::endl;
        // Receive
        if ((reducerRank % indexPivot) == 0){
            sender = reducerRank + previousPivot;
            // Does the sender exist?
            if (sender < numReducers){
                messageBuffer = new char[bufferSize];
                std::cout << mpiRank << " 7.1.1 - " << reducerRank << " receiving: " << (sender + reducerStartRank) << std::endl;
                MPI_Recv(messageBuffer, bufferSize, MPI_CHAR, (sender + reducerStartRank), REDUCE_MESSAGE, MPI_COMM_WORLD, &status);
                std::cout << mpiRank << " 7.2 - " << hostname << " - received from: " << (sender + reducerStartRank) << std::endl;
                extractPointPairsFromMessage(messageBuffer, receivedSortedPairs);
                delete[] messageBuffer;
                mergeSortedArrays(receivedSortedPairs);
            }
        // Send
        } else if ((reducerRank % previousPivot) == 0 || previousPivot == 1){
            receiver = reducerRank - previousPivot;
            message = serializePairList(sortedPairs);
            messageBuffer = new char[bufferSize];
            strcpy(messageBuffer, message->c_str());
            std::cout << mpiRank << " 7.3 - " << hostname << " - sending to: " << (receiver + reducerStartRank) << std::endl;
            int ierr = MPI_Send(messageBuffer, bufferSize, MPI_CHAR, (receiver + reducerStartRank), REDUCE_MESSAGE, MPI_COMM_WORLD);
            this->checkForSendErrors(ierr, "uff::TopKSimilarityMapRedApp::runFinalBinaryReduce");
            delete[] messageBuffer;
        } else {
            break;
        }
        previousPivot = indexPivot;
        indexPivot = indexPivot * 2;
    }
}

std::string* uff::TopKSimilarityMapRedApp::serializePointList(std::vector<uff::Point*>* pointList){
    std::vector<Point*>::iterator it;
    std::string* message = new std::string();
    std::string separator = "";
    for (it = pointList->begin(); it < pointList->end(); it++) {
        message->append(separator);
        message->append((*it)->serialize());
        separator="|";
    }
    return message;
}

std::string* uff::TopKSimilarityMapRedApp::serializePairList(std::vector<uff::PointPair*>* pairList){
    std::vector<PointPair*>::iterator it;
    std::string* message = new std::string();
    std::string separator = "";
    for (it = pairList->begin(); it < pairList->end(); it++) {
        message->append(separator);
        message->append((*it)->serialize());
        separator="|";
    }
    return message;
}

void uff::TopKSimilarityMapRedApp::mergeSortedArrays(std::vector<uff::PointPair*>* receivedSortedPairs){
    std::vector<PointPair*>* localPairs = new std::vector<PointPair*>();
    std::vector<PointPair*>::iterator itLocal = sortedPairs->begin();
    std::vector<PointPair*>::iterator itRemote = receivedSortedPairs->begin();
    while(itLocal < sortedPairs->end() || itRemote < receivedSortedPairs->end()){
        if (itRemote >= receivedSortedPairs->end()){
            localPairs->push_back((* itLocal));
            itLocal++;
        } else if (itLocal >= sortedPairs->end()){
            localPairs->push_back( (* itRemote ));
            itRemote++;
        } else if ((*itRemote)->similarity < (*itLocal)->similarity) {
            localPairs->push_back((*itRemote));
            itRemote++;
        } else{
            localPairs->push_back((*itLocal));
            itLocal++;
        }
    }
    for (int i=localPairs->size() - 1; i>=this->k; i--){
        delete(localPairs->at(i));
        localPairs->pop_back();
    }

    delete(sortedPairs);
    sortedPairs = localPairs;
}

void uff::TopKSimilarityMapRedApp::writeResultToFile(){

    std::ofstream resultFile;
    std::string outputFileName = this->outpuFolder;
    outputFileName.append("/output");
    resultFile.open(outputFileName.c_str());
    if (!resultFile.is_open()){
        std::cout << "Couldn't open output file " << outputFileName << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
        exit(2);
    }
    for (int i=0; i < sortedPairs->size(); i++) {
        resultFile << i << "\t" << sortedPairs->at(i)->serialize() << std::endl;
    }
    resultFile.close();
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Result wrote  @" << this->hostname << ":" <<  outputFileName << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
}

void uff::TopKSimilarityMapRedApp::receiveFinalResultFromReduce(){
    char *messageBuffer = new char[bufferSize];
    sortedPairs = new std::vector<uff::PointPair*>();
    MPI_Status status;
    std::cout << mpiRank << " 8.2 - Waiting for final message message from " << reducerStartRank << std::endl;
    MPI_Recv(messageBuffer, bufferSize, MPI_CHAR, reducerStartRank, REDUCE_MESSAGE, MPI_COMM_WORLD, &status);
    extractPointPairsFromMessage(messageBuffer, sortedPairs);
    delete[] messageBuffer;
}

void uff::TopKSimilarityMapRedApp::sendResultToMaster(){
    std::string* message = serializePairList(sortedPairs);
    char *messageBuffer = new char[bufferSize];
    strcpy(messageBuffer, message->c_str());
    std::cout << mpiRank << " 8.1 - sending final result to master" << std::endl;
    int ierr = MPI_Send(messageBuffer, bufferSize, MPI_CHAR, 0, REDUCE_MESSAGE, MPI_COMM_WORLD);
    this->checkForSendErrors(ierr, "uff::TopKSimilarityMapRedApp::sendResultToMaster");
    delete[] messageBuffer;

}

// MAP --------------------------------------------------------
void uff::TopKSimilarityMapRedApp::pointToPairMapper(long index, Point* pt) {

    long pid = (index % (m)) + 1; // [1, m]

    pt->setIndex(index);
    pt->setPid(pid);

    int reduceIndex = 0;

    std::string message;

    std::stringstream key;
    for (long i = 1; i <= m; i++) {
        key.str(std::string());
        // index(i,j) = ((i-1) * m) + j - sum(i-1) - 1
        if (i <= pid) {
            reduceIndex = ((i - 1) * m) + pid - ((i * (i - 1)) / 2) - 1;
            key << "[" << i << ", " << pid << "]";
        } else {
            reduceIndex = ((pid - 1) * m) + i - ((pid * (pid - 1)) / 2) - 1;
            key << "[" << pid << ", " << i << "]";
        }
//        std::cout << mpiRank << " - " << key.str() << "Index: " << index << " pid: " << pid << std::endl;
        registerMessageToReduce(reduceIndex, key.str(), pt->serialize());
    }
}

// Reduce ------------------------------------------------
void uff::TopKSimilarityMapRedApp::topKDistanceReducer(uff::KeyPair* key,
        std::vector<uff::Point*> *points) {

    std::vector<Point*> pointList_i;
    std::vector<Point*> pointList_j;

    sortedPairs = new std::vector<PointPair*>();

    try{
        // Set the cartesian product between points with Pid1 and points with Pid2
        for (std::vector<int>::size_type i = 0; i < points->size(); i++) {
            if (points->at(i) == NULL){
                std::cerr << mpiRank << " ### ERROR: " << i << std::endl;
                continue;
            }
    //        std::cout << mpiRank << " - " << pointList[i]->getPid() << std::endl;
            if (pointList[i]->getPid() == key->value1) {
                pointList_i.push_back(points->at(i));
            }
            if (pointList[i]->getPid() == key->value2) {
                pointList_j.push_back(points->at(i));
            }
        }
    } catch (std::bad_alloc& ba) {
        std::cerr << mpiRank << " topKDistanceReducer creating vector bad_alloc caught: " << ba.what() << std::endl;
        return;
    }
//    std::cout << mpiRank << " - " << key->serialize() << " i size: " << serializePointList(&pointList_i)->c_str() << std::endl;
//    std::cout << mpiRank << " - " << key->serialize() << " j size: " << serializePointList(&pointList_j)->c_str() << std::endl;

    std::vector<PointPair*>::iterator it;

    PointPair* pair;
    std::vector<int>::size_type i, j;
    int startJ=0;
    bool inserted;
    try{
        for (i=0; i<pointList_i.size(); i++){
            startJ=0;
            if (key->value1==key->value2) startJ=i + 1;
            for (j = startJ; j < pointList_j.size();    j++) {
                if (pointList_i.at(i)->equals(pointList_j.at(j))){
                    continue;
                }
    //            std::cout << mpiRank << " similarity " << pointList_i[i]->serialize() << " - " <<   pointList_j[j]->serialize() << std::endl;
                double similarity = uff::Calc::calculateSimilarity(pointList_i[i], pointList_j[j]);
                pair = new uff::PointPair(pointList_i[i], pointList_j[j], similarity);
    //            std::cout << ".";
                // Add pair to the fixed size sorted list
                inserted = false;
                if (sortedPairs->size() == 0) {
                    sortedPairs->push_back(pair);
                    inserted = true;
                } else {
                    for (it = sortedPairs->begin(); it < sortedPairs->end(); it++) {
                        if (pair->similarity < (*it)->similarity) {
                            sortedPairs->insert(it, pair);
                            inserted = true;
                            break;
                        } else if (pair->similarity == (*it)->similarity
                                && pair->pointA->equals((*it)->pointA)
                                && pair->pointB->equals((*it)->pointB)){
                            break;
                        }
                    }
                }
                if (!inserted){
                    delete(pair);
                }
                if (sortedPairs->size() > k) {
                    pair = sortedPairs->back();
                    sortedPairs->pop_back();
                    delete(pair);
                }
            }
        }
    } catch (std::bad_alloc& ba) {
        std::cerr << mpiRank << " topKDistanceReducer running bad_alloc caught: " << ba.what() << std::endl;
        return;
    }
}

void uff::TopKSimilarityMapRedApp::checkForSendErrors(int ierr, std::string msg){
    int errclass;
    char err_buffer[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_class(ierr, &errclass);
    if (errclass == MPI_ERR_RANK) {
        std::cout << mpiRank << " error sending message: " << msg << std::endl;
        MPI_Error_string(ierr, err_buffer, &resultlen);
        std::cerr << err_buffer << std::endl;
        MPI_Finalize(); /* abort*/
    }
}

void uff::TopKSimilarityMapRedApp::log(std::string msg){
    std::cout << this->mpiRank << " " << msg << std::endl;
}

int main(int argc, char *argv[]) {
    int mpiSize;
    int mpiRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);

    // Read input params
    if (argc <= 5) {
        if (mpiRank == 0) {
            std::cerr << "Usage: " << argv[0] << " <k> <m> <n> <inputFolder> <outputFolder>" << std::endl;
            std::cerr << "  k - top k pairs" << std::endl;
            std::cerr << "  m - numReducers = (m * (m + 1))/2" << std::endl;
            std::cerr << "  n - number of files (same as number of mappers)" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    int k = atoi(argv[1]);
    int m = atoi(argv[2]);
    int numFiles = atoi(argv[3]);
    std::string inputFolder = argv[4];
    std::string outputFolder = argv[5];

    uff::TopKSimilarityMapRedApp* app = new uff::TopKSimilarityMapRedApp(mpiSize,
            mpiRank);

    app->run(k, m, numFiles, inputFolder, outputFolder);
    app->log(" 9 - ##################################### finalizing...");
    MPI_Finalize();
    return 0;
}
