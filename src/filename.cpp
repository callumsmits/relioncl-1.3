/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/filename.h"
#include "src/funcs.h"

// Constructor with root, number and extension .............................
void FileName::compose(const std::string &str, long int no, const std::string &ext, int numberlength)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[numberlength+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(numberlength)+
                         (std::string)"d";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this += aux_str;
    }

    if (ext != "")
        *this += (std::string)"." + ext;
}

// Constructor: prefix number and filename, mainly for selfiles..
void FileName::compose(long int no , const std::string &str, int numberlength)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[numberlength+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(numberlength)+
                         (std::string)"d@";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this = aux_str + str;
    }
    else
        *this = str;


}

// Is in stack ............................................................
bool FileName::isInStack() const
{
    return find("@") != std::string::npos;
}

// Decompose ..............................................................
void FileName::decompose(long int &no, std::string &str) const
{
    size_t idx = find('@');
    if(idx != std::string::npos)
    {
        no = textToInteger(substr(0,idx));
        str = substr(idx+1,length()-idx);
    }
    else{
      no=-1;
      str = *this;
    }
}

// Convert to lower case characters .........................................
FileName FileName::toLowercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = tolower(result[i]);
    return result;
}

// Convert to upper case characters .........................................
FileName FileName::toUppercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = toupper(result[i]);
    return result;
}

// Is substring present?
bool FileName::contains(const std::string& str) const
{
    int point = rfind(str);
    if (point > -1)
        return true;
    else
        return false;
}

// Get substring before first instance of str
FileName FileName::beforeFirstOf(const std::string& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring before last instance of str
FileName FileName::beforeLastOf(const std::string& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring after first instance of str
FileName FileName::afterFirstOf(const std::string& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get substring after last instance of str
FileName FileName::afterLastOf(const std::string& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get the base name of a filename .........................................
std::string FileName::getBaseName() const
{
    std::string basename = "";
    std::string myname = *this;
    int myindex = 0;
    for (int p = myname.size() - 1; p >= 0; p--)
    {
        if (myname[p] == '/')
        {
            myindex = p + 1;
            break;
        }
    }
    for (int p = myindex; p < myname.size(); p++)
    {
        if (myname[p] != '.')
            basename += myname[p];
        else
            break;
    }
    return basename;
}

// Get the extension of a filename .........................................
std::string FileName::getExtension() const
{
    int last_point = find_last_of(".");
    if (last_point == -1)
        return "";
    else
        return substr(last_point + 1);
}

// Init random .............................................................
void FileName::initRandom(int length)
{
    randomize_random_generator();
    *this = "";
    for (int i = 0; i < length; i++)
        *this += 'a' + FLOOR(rnd_unif(0, 26));
}

// Add at beginning ........................................................
FileName FileName::addPrefix(const std::string &prefix) const
{
    FileName retval = *this;
    int skip_directories = find_last_of("/") + 1;
    return retval.insert(skip_directories, prefix);
}

// Add at the end ..........................................................
FileName FileName::addExtension(const std::string &ext) const
{
    if (ext == "")
        return *this;
    else
    {
        FileName retval = *this;
        retval = retval.append((std::string)"." + ext);
        return retval;
    }
}

// Remove last extension ...................................................
FileName FileName::withoutExtension() const
{
    FileName retval = *this;
    return retval.substr(0, rfind("."));
}

// Insert before extension .................................................
FileName FileName::insertBeforeExtension(const std::string &str) const
{
    int point = -1;
    bool done = false;
    do
    {
        point = find(".", point + 1);
        if (point == -1)
        {
            point = length();
            done = true;
        }
        else if (point == length() - 1)
            done = true;
        else if ((*this)[point+1] == '.' || (*this)[point+1] == '/')
            done = false;
        else
            done = true;
    }
    while (!done);
    FileName retval = *this;
    return retval.insert(point, str);
}

// Remove an extension wherever it is ......................................
FileName FileName::removeExtension(const std::string &ext) const
{
    int first = find((std::string)"." + ext);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(first, 1 + ext.length());
    }
}

// Remove all extensions....................................................
FileName FileName::removeAllExtensions() const
{
    int first = rfind("/");
    first = find(".", first + 1);
    if (first == -1)
        return *this;
    else
        return substr(0, first);
}

FileName FileName::getFileFormat() const
{
    int first;
    FileName result;
    if (find("#") != std::string::npos)
        return "raw";
    else if ( (first = rfind(":"))!=std::string::npos)
        result = substr(first + 1) ;
    else if ( (first = rfind("."))!=std::string::npos)
        result = substr(first + 1);
    else
        result="spi";
    return result.toLowercase();

}

FileName FileName::removeFileFormat() const
{
    if ( find("#", 0) > -1 )
        REPORT_ERROR("Not implemented for raw data");
    size_t found=rfind(":");
    if (found!=std::string::npos)
        return substr(0, found);
    return *this;
}

bool FileName::isStarFile() const
{
    //file names containing @, : or % are not metadatas
    size_t found=this->find('@');
    if (found!=std::string::npos)
        return false;
    found=this->find(':');
    if (found!=std::string::npos)
        return false;
    found=this->find('#');
    if (found!=std::string::npos)
        return false;

    FileName ext = getFileFormat();
    if (ext=="star")
    {
        return true;
    }
    else
    {
    	return false;
    }
}

// Substitute one extension by other .......................................
FileName FileName::substituteExtension(const std::string &ext1,
                                       const std::string &ext2) const
{
    int first = find((std::string)"." + ext1);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.replace(first, 1 + ext1.length(), (std::string)"." + ext2);
    }
}

// Remove a substring ......................................................
FileName FileName::without(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(pos, str.length());
    }
}

// Remove until prefix .....................................................
FileName FileName::removeUntilPrefix(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(0, pos + str.length());
    }
}

// Remove directories ......................................................
FileName FileName::removeDirectories(int keep) const
{
    int last_slash = rfind("/");
    int tokeep = keep;
    while (tokeep > 0)
    {
        last_slash = rfind("/", last_slash - 1);
        tokeep--;
    }
    if (last_slash == -1)
        return *this;
    else
        return substr(last_slash + 1, length() - last_slash);
}
void FileName::copyFile(const FileName & target) const
{
    std::ifstream f1 (this->c_str(), std::fstream::binary);
    std::ofstream f2 (target.c_str(),std::fstream::trunc|std::fstream::binary);
    f2<<f1.rdbuf();
}

int FileName::globFiles(std::vector<FileName> &files) const
{
	glob_t glob_result;
	glob((*this).c_str(), GLOB_TILDE, NULL, &glob_result);
	files.clear();
	for(unsigned  int  i = 0; i < glob_result.gl_pathc; ++i)
	{
		files.push_back(std::string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);
	return files.size();
}

bool exists(const FileName &fn)
{
    FILE *aux;
    if ((aux = fopen(fn.c_str(), "r")) == NULL)
        return false;
    fclose(aux);
    return true;
}
