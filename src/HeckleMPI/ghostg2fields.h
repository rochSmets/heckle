#ifndef GHOSTCOMMSFIELDS_H
#define GHOSTCOMMSFIELDS_H


#include <misc.h>
//#include <ghosts.h>



#define GHOSTFIELD_E 0
#define GHOSTFIELD_F 1
#define GHOSTFIELD_J 2


typedef struct ghostg2fields_s GhostG2Fields;





GhostG2Fields* GhostG2FieldsInit(const STX * const sx, int fieldID);



void GhostG2FieldsDelete(GhostG2Fields* self);



void GhostG2FieldsSendRecv(GhostG2Fields *self, const STX* const sx, ST2 * s2);




#endif // GHOSTCOMMSFIELDS_H

