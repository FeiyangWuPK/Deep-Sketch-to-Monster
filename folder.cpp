#include "folder.h"
#include <QFileDialog>

std::string lf = "";

const std::string & localfolder(){
    if (!lf.length()){
        std::string localfolder_linux = QDir::currentPath().toStdString();
        for (int i = 0; i < localfolder_linux.length(); ++i){
            if (localfolder_linux[i] == '/'){
                lf += '\\';
            } else lf += localfolder_linux[i];
        }
    }
    return lf;
}



