#pragma once

class PatternID
{
public:
	PatternID(short ,short);
	PatternID(PatternID* ,short );
	~PatternID();
	char* toString(void) const;
	short * getValue() const;
    short getSize() const;
	void print()const;
	short getIDOfNodeAt(short index) const;
	bool containsNode(short) const;
	short getNumberOfOverlap(PatternID* otherPatternID) const;
	static int counter;
private:
	short * value;
	short size;
};
