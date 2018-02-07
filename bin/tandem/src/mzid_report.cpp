/*
 Copyright (C) 2014 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2014-09-25
/*
 * the mreport object provides the functionality to export all of the information collected
 * in mprocess to an XML file. The file path was defined in the original XML input parameter
 * file. The file produced by mreport can be used as an input parameter file to repeat
 * a protein modeling session at a later time.
 */

#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
//#include <regex>
#include "saxsaphandler.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"

// create the mzid_report object and load the UniMod definitions

mzid_report::mzid_report(mscore& score)
	: m_Score(score)
{
	load_unimod();
}

mzid_report::~mzid_report(void)
{
}

// Open the ofstream object using the file name supplied in the input file and
// add the header portion of the mzIndentML file
// Should be called initially in an use of the object

bool mzid_report::start(XmlParameter &_x,XmlParameter &_p)
{
	string strKey = "output, path";
	string strValue;
	_x.get(strKey,strValue);
	if(strValue.size() == 0)
		return false;
	m_strPath = strValue;
	m_strPath.append(".mzid");
	m_ofOut.open(m_strPath.c_str());
	if(m_ofOut.fail())	{
		return false;
	}
	strKey = "process, version";
	_p.get(strKey,strValue);
	string strVersion = strValue;
	strKey = "process, start time";
	_p.get(strKey,strValue);
	string strTime = strValue;
	strKey = "output, mzid decoy DB accession regexp";
	_x.get(strKey,strValue);
	if(strValue.length() == 0)	{
		m_strRegex = "^XXX";
	}
	else	{
		m_strRegex = strValue;
	}
	return add_header(strVersion,strTime);
}

// Finish the XML output and close the ofstream object
// Should be called last in any use of the object

bool mzid_report::end(void)
{
	if(m_ofOut.fail())	{
		return false;
	}
	m_ofOut << "</MzIdentML>\n";
	m_ofOut.close();
	return true;
}

// Output the mzIndentML header information

bool mzid_report::add_header(string &_s,string &_t)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}

	m_ofOut << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	m_ofOut << "<MzIdentML id=\"X! Tandem\" version=\"1.1.0\" xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psidev.info/psi/pi/mzIdentML/1.1 http://www.psidev.info/files/mzIdentML1.1.0.xsd\" creationDate=\"" << _t << "\" >\n";
	m_ofOut << "<cvList xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">\n";
	m_ofOut << "    <cv id=\"PSI-MS\" uri=\"http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.30.0\" fullName=\"PSI-MS\"/>\n";
	m_ofOut << "    <cv id=\"UNIMOD\" uri=\"http://www.unimod.org/obo/unimod.obo\" fullName=\"UNIMOD\"/>\n";
	m_ofOut << "    <cv id=\"UO\" uri=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\" fullName=\"UNIT-ONTOLOGY\"/>\n";
	m_ofOut << "</cvList>\n";
	m_ofOut << "<AnalysisSoftwareList xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">\n";
	m_ofOut << "    <AnalysisSoftware version=\"" << _s << "\" name=\"X! Tandem\" id=\"ID_software\">\n";
	m_ofOut << "        <SoftwareName>\n";
	m_ofOut << "            <cvParam accession=\"MS:1001476\" cvRef=\"PSI-MS\" name=\"X\\! Tandem\"/>\n";
	m_ofOut << "        </SoftwareName>\n";
	m_ofOut << "    </AnalysisSoftware>\n";
	m_ofOut << "</AnalysisSoftwareList>\n";
	return true;
}

// Create the SequenceCollection portion of the mzIndentML file

bool mzid_report::sequence_collection(vector<mspectrum> &_vs,vector<string> &_vp)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	m_ofOut << "<SequenceCollection xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">\n";
	add_dbsequence(_vs);
	add_peptides(_vs);
	add_peptide_evidence(_vs);
	m_ofOut << "</SequenceCollection>\n";
	return true;
}

// Create the AnalysisCollection portion of the mzIdentML file

bool mzid_report::analysis_collection(vector<mspectrum> &_vs,vector<string> &_vp)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	m_ofOut << "<AnalysisCollection xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">";
	m_ofOut << "\t<SpectrumIdentification spectrumIdentificationList_ref=\"SI_LIST_1\" spectrumIdentificationProtocol_ref=\"SearchProtocol_1\" id=\"SpecIdent_1\">\n";
	m_ofOut << "\t\t<InputSpectra spectraData_ref=\"SID_1\"/>\n";
	size_t a = 0;
	while(a < _vp.size())	{
		m_ofOut << "\t\t<SearchDatabaseRef searchDatabase_ref=\"SearchDB_" << a <<"\"/>\n";
		a++;
	}
	m_ofOut << "\t</SpectrumIdentification>\n";
	m_ofOut << "</AnalysisCollection>\n";
	return true;
}

// Create the AnalysisProtocolCollection portion of the mzIdentML file

bool mzid_report::analysis_protocol_collection(XmlParameter &_x)
{
	m_ofOut << "<AnalysisProtocolCollection xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">\n";
	m_ofOut << "<SpectrumIdentificationProtocol analysisSoftware_ref=\"ID_software\" id=\"SearchProtocol_1\">\n";
	m_ofOut << "\t<SearchType>\n";
    m_ofOut << "\t\t<cvParam accession=\"MS:1001083\" cvRef=\"PSI-MS\" name=\"ms-ms search\"/>\n";
    m_ofOut << "\t</SearchType>\n";
    m_ofOut << "\t<ModificationParams>\n";

	string strKey;
	string strValue;
	strKey = "residue, modification mass";
	_x.get(strKey,strValue);
	analysis_mods(strValue,true);
	strKey = "residue, potential modification mass";
	_x.get(strKey,strValue);
	analysis_mods(strValue,false);
	strKey = "refine, potential modification mass";
	_x.get(strKey,strValue);
	analysis_mods(strValue,false);
	strKey = "refine, potential modification mass 1";
	_x.get(strKey,strValue);
	analysis_mods(strValue,false);
	
	m_ofOut << "\t</ModificationParams>\n";
    m_ofOut << "\t<Enzymes>\n";
    m_ofOut << "\t<Enzyme missedCleavages=\"1000\" semiSpecific=\"false\" id=\"Tryp\">\n";
    m_ofOut << "\t<EnzymeName>\n";
    m_ofOut << "\t\t<cvParam accession=\"MS:1001251\" cvRef=\"PSI-MS\" name=\"Trypsin\"/>\n";
    m_ofOut << "\t</EnzymeName>\n";
    m_ofOut << "\t</Enzyme>\n";
    m_ofOut << "\t</Enzymes>\n";
    m_ofOut << "\t<ParentTolerance>\n";
    m_ofOut << "\t<cvParam accession=\"MS:1001412\" cvRef=\"PSI-MS\" unitCvRef=\"UO\" unitName=\"parts per million\" unitAccession=\"UO:0000169\" value=\"20.0\" name=\"search tolerance plus value\"/>\n";
    m_ofOut << "\t<cvParam accession=\"MS:1001413\" cvRef=\"PSI-MS\" unitCvRef=\"UO\" unitName=\"parts per million\" unitAccession=\"UO:0000169\" value=\"20.0\" name=\"search tolerance minus value\"/>\n";
    m_ofOut << "\t</ParentTolerance>\n";
    m_ofOut << "\t<Threshold>\n";
    m_ofOut << "\t<cvParam accession=\"MS:1001494\" cvRef=\"PSI-MS\" name=\"no threshold\"/>\n";
    m_ofOut << "\t</Threshold>\n";
    m_ofOut << "\t</SpectrumIdentificationProtocol>\n";
    m_ofOut << "</AnalysisProtocolCollection>\n";
	return true;
}

// Create the DataCollection portion of the mzIndentML file

bool mzid_report::data_collection(vector<mspectrum> &_vs,vector<string> &_vp,XmlParameter &_x)
{
	m_ofOut << "<DataCollection xmlns=\"http://psidev.info/psi/pi/mzIdentML/1.1\">\n";

	add_inputs(_vs,_vp,_x);
	add_analysis(_vs,_vp,_x);

	m_ofOut << "</DataCollection>\n";

	return true;
}

// Create the AnalysisData portion of the mzIndentML file
// Used by the data_collection method

bool mzid_report::add_analysis(vector<mspectrum> &_vs,vector<string> &_vp,XmlParameter &_x)
{
	m_ofOut << "<AnalysisData>\n";
    m_ofOut << "<SpectrumIdentificationList id=\"SI_LIST_1\">\n";
    m_ofOut << " <FragmentationTable>\n";
    m_ofOut << "<Measure id=\"Measure_MZ\">\n";
    m_ofOut << "<cvParam accession=\"MS:1001225\" cvRef=\"PSI-MS\" unitCvRef=\"PSI-MS\" unitName=\"m/z\" unitAccession=\"MS:1000040\" name=\"product ion m/z\"/>\n";
    m_ofOut << "</Measure>\n";
	m_ofOut << "</FragmentationTable>\n";
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	char *pLine = new char[1024];
	int iScan = 0;
	string strScan;
	size_t tPos = 0;
	double	dProton = 1.007276;
	vector<string> vstrPeptides;
	vector<size_t> vtB;
	vector<size_t> vtC;
	string strPep;
	while(a < _vs.size())	{
		if(_vs[a].m_vseqBest.empty())	{
			a++;
			continue;
		}
		iScan = (int)(_vs[a].m_tId);
		tPos = _vs[a].m_strDescription.find("scan=");
		if(tPos != _vs[a].m_strDescription.npos)	{
			tPos += 5;
			strScan = _vs[a].m_strDescription.substr(tPos,_vs[a].m_strDescription.length() - tPos + 1);
			iScan = atoi(strScan.c_str());
		}
        m_ofOut << "<SpectrumIdentificationResult spectraData_ref=\"SID_1\" ";
		m_ofOut << "spectrumID=\"" << _vs[a].m_strDescription <<"\" ";
		m_ofOut << "id=\"SIR_" << iScan << "\">\n";
 		b = 0;
		vstrPeptides.clear();
		vtC.clear();
		vtB.clear();
		while(b < _vs[a].m_vseqBest.size())	{
			c = 0;
			while(c < _vs[a].m_vseqBest[b].m_vDomains.size())	{
				sprintf(pLine,"%i_%i_%i", (int)_vs[a].m_tId,(int)(b+1),(int)(c+1));
				if(m_mapIdsToPeptides.find(pLine) == m_mapIdsToPeptides.end())	{
					c++;
					continue;
				}
				strPep = m_mapIdsToPeptides.find(pLine)->second;
				vstrPeptides.push_back(strPep);
				vtB.push_back(b);
				vtC.push_back(c);
				c++;
			}
			b++;
		}
		size_t tOut = 0;
		string strCurrentPeptide;
		size_t x = 0;
		string strPeptide;
		while(tOut < vstrPeptides.size())	{
			x = 0;
			while(x < vstrPeptides.size() && vstrPeptides[x].empty())	{
				x++;
			}
			if(x == vstrPeptides.size())	{
				break;
			}
			strCurrentPeptide = vstrPeptides[x];
			x = 0;
			b = 0;
			m_ofOut << "<SpectrumIdentificationItem passThreshold=\"true\" ";
			m_ofOut << "rank=\"1\" peptide_ref=\"Pep_" << strCurrentPeptide << "\" ";
			sprintf(pLine,"%.4lf",dProton + (_vs[a].m_vseqBest[vtB[x]].m_vDomains[vtC[x]].m_dMH - dProton)/(double)(_vs[a].m_fZ));
			m_ofOut << "calculatedMassToCharge=\"" << pLine << "\" ";
			sprintf(pLine,"%.4lf",dProton + (_vs[a].m_dMH - dProton)/(double)(_vs[a].m_fZ));
			m_ofOut << "experimentalMassToCharge=\"" << pLine << "\" ";
			m_ofOut << "chargeState=\"" << (int)_vs[a].m_fZ << "\" ";
			m_ofOut << "id=\"SII_" << _vs[a].m_tId << "_" << x << "\">\n";
			while(b < _vs[a].m_vseqBest.size())	{
				c = 0;
				while(c < _vs[a].m_vseqBest[b].m_vDomains.size())	{
					sprintf(pLine,"%i_%i_%i",(int)_vs[a].m_tId,(int)(b+1),(int)(c+1));
					strPeptide = pLine;
					if(strCurrentPeptide == vstrPeptides[x])	{
						m_ofOut << "<PeptideEvidenceRef peptideEvidence_ref=\"PepEv_" << strPeptide << "\"/>\n";
						tOut++;
						vstrPeptides[x].clear();
					}
					x++;
					c++;
				}
				b++;
			}
			sprintf(pLine,"%.1e",_vs[a].m_hHyper.expect(m_Score.hconvert(_vs[a].m_vseqBest[0].m_vDomains[0].m_fHyper)));
			m_ofOut << "<cvParam accession=\"MS:1001330\" cvRef=\"PSI-MS\" value=\"" << pLine << "\" name=\"X\\!Tandem:expect\"/>\n";
			m_Score.report_score(pLine, _vs[a].m_vseqBest[0].m_vDomains[0].m_fHyper);
			m_ofOut << "<cvParam accession=\"MS:1001331\" cvRef=\"PSI-MS\" value=\"" << pLine << "\" name=\"X\\!Tandem:hyperscore\"/>\n";
			m_ofOut << "</SpectrumIdentificationItem>\n";
		}
		m_ofOut << "<cvParam accession=\"MS:1001115\" cvRef=\"PSI-MS\" value=\"" << iScan << "\" name=\"scan number(s)\"/>\n";
		m_ofOut << "</SpectrumIdentificationResult>\n";

		a++;
	}
	m_ofOut << "</SpectrumIdentificationList>\n</AnalysisData>\n";
	delete pLine;
	return true;
}

// Create the Inputs portion of the mzIndentML file
// Used by the data_collection method

bool mzid_report::add_inputs(vector<mspectrum> &_vs,vector<string> &_vp,XmlParameter &_x)
{
	m_ofOut << "<Inputs>\n";
	string strKey = "spectrum, path";
	size_t a = 0;
	string strValue;
	size_t tPos = 0;
	while(a < _vp.size())	{
        m_ofOut << "<SearchDatabase location=\"" << _vp[a] << "\" id=\"SearchDB_" << a << "\">\n";
        m_ofOut << "\t<FileFormat>\n";
        m_ofOut << "\t\t<cvParam accession=\"MS:1001348\" cvRef=\"PSI-MS\" name=\"FASTA format\"/>\n";
        m_ofOut << "\t</FileFormat>\n";
        m_ofOut << "\t<DatabaseName>\n";
		tPos = _vp[a].find_last_of('/');
		strValue = _vp[a];
		if(tPos != _vp[a].npos)	{
				tPos++;
				strValue = _vp[a].substr(tPos,_vp[a].size() - tPos + 1);
		}
		else	{
			tPos = _vp[a].find_last_of('\\');
			if(tPos != _vp[a].npos)	{
					tPos++;
					strValue = _vp[a].substr(tPos,_vp[a].size() - tPos + 1);
			}
		}

        m_ofOut << "\t\t<userParam name=\"" << strValue << "\"/>\n";
        m_ofOut << "\t</DatabaseName>\n";
        m_ofOut << "\t<cvParam accession=\"MS:1001197\" cvRef=\"PSI-MS\" name=\"DB composition target+decoy\"/>\n";
		strKey = "output, mzid decoy DB accession regexp";
		_x.get(strKey,strValue);
		if(strValue.size() > 0)	{
	        m_ofOut << "\t<cvParam accession=\"MS:1001283\" cvRef=\"PSI-MS\" value=\"" << strValue << "\" name=\"decoy DB accession regexp\"/>\n";
		}
		else {
	        m_ofOut << "\t<cvParam accession=\"MS:1001283\" cvRef=\"PSI-MS\" value=\"^XXX\" name=\"decoy DB accession regexp\"/>\n";
		}
        m_ofOut << "\t<cvParam accession=\"MS:1001195\" cvRef=\"PSI-MS\" name=\"decoy DB type reverse\"/>\n";
        m_ofOut << "\t</SearchDatabase>\n";
		a++;
	}
	strKey = "spectrum, path";
	_x.get(strKey,strValue);
    m_ofOut << "\t<SpectraData location=\"" << strValue << "\" ";
	tPos = strValue.rfind("/");
	if(tPos != strValue.npos)	{
			tPos++;
			strValue = strValue.substr(tPos,strValue.size() - tPos + 1);
	}
	else	{
		tPos = strValue.rfind("\\");
		if(tPos != strValue.npos)	{
				tPos++;
				strValue = strValue.substr(tPos,strValue.size() - tPos + 1);
		}
	}
    m_ofOut << "name=\"" << strValue << "\" id=\"SID_1\">\n";
    m_ofOut << "\t<FileFormat>\n";
    m_ofOut << "\t\t<cvParam accession=\"MS:1000058\" cvRef=\"PSI-MS\" name=\"mzML file\"/>\n";
    m_ofOut << "\t</FileFormat>\n";
    m_ofOut << "\t<SpectrumIDFormat>\n";
    m_ofOut << " \t\t<cvParam accession=\"MS:1000768\" cvRef=\"PSI-MS\" name=\"Thermo nativeID format\"/>\n";
    m_ofOut << "\t</SpectrumIDFormat>\n";
     m_ofOut << "\t</SpectraData>\n";

	 m_ofOut << "</Inputs>\n";

	return true;
}

// Create the residue modification portion of the mzIndentML file
// Used by the analysis_protocol_collection method

bool mzid_report::analysis_mods(string &_v,bool _b)
{
	string strValue = _v;
	char *pMass = new char[256];
	char cRes = '\0';
	size_t tA = 0;
	size_t tB = 0;
	size_t tLength = strValue.length();
	string strUni;
	string strUniDescription;
	while(tA < tLength)	{
		strcpy(pMass,"");
		tB = 0;
		while(strValue[tA] != '@' && tA < tLength)	{
			pMass[tB] = strValue[tA];
			tB++;
			tA++;
		}
		pMass[tB] = '\0';
		tA++;
		cRes = strValue[tA];
        m_ofOut << "\t\t<SearchModification residues=\"" << cRes << "\" massDelta=\"" << atof(pMass) << "\" fixedMod=\"";
		if(_b)	{
			m_ofOut << "true\">\n";
		}
		else	{
			m_ofOut << "false\">\n";
		}
		if(get_unimod(atof(pMass),strUni,strUniDescription))	{
	        m_ofOut << "\t\t\t<cvParam accession=\"" << strUni << "\" cvRef=\"UNIMOD\" name=\"" << strUniDescription << "\"/>\n";
		}
        m_ofOut << "\t\t</SearchModification>\n";
		while(strValue[tA] != ',' && tA < tLength)	{
			tA++;
		}
		tA++;
	}
	delete pMass;
	return true;
}

// Create the DBSequence entries in the mzIndentML file
// Used by the sequence_collection method

bool mzid_report::add_dbsequence(vector<mspectrum> &_vs)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	set<size_t> setUids;
	set<size_t>::iterator itUid;
	const size_t tSize = _vs.size();
	size_t a = 0;
	size_t b = 0;
	size_t tUid = 0;
	string strLabel;
	string strDesc;
	while(a < tSize)	{
		if(_vs[a].m_vseqBest.empty())	{
			a++;
			continue;
		}
		b = 0;
		while(b < _vs[a].m_vseqBest.size())	{
			tUid = _vs[a].m_vseqBest[b].m_tUid;
			itUid = setUids.find(tUid);
			if(itUid == setUids.end())	{
				parse_description(_vs[a].m_vseqBest[b].m_strDes,strLabel,strDesc);
				m_ofOut << "<DBSequence accession=\"" << strLabel.c_str() << "\" searchDatabase_ref=\"SearchDB_" << _vs[a].m_vseqBest[b].m_siPath << "\" length=\"" << (unsigned long)_vs[a].m_vseqBest[b].m_strSeq.size() << "\" id=\"DBSeq" << tUid << "\">\n";
				m_ofOut << "<cvParam accession=\"MS:1001088\" cvRef=\"PSI-MS\" value=\"" << strDesc.c_str() << "\" name=\"protein description\"/>\n";
				m_ofOut << "</DBSequence>\n";
				setUids.insert(tUid);
			}
			b++;
		}
		a++;
	}
	return true;
}

// Create the PeptideEvidence entries in the mzIndentML file
// Used by the sequence_collection method

bool mzid_report::add_peptide_evidence(vector<mspectrum> &_vs)
{
	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	const size_t tSize = _vs.size();
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	size_t tUid = 0;
	string strPre;
	string strPost;
	size_t tStart = 0;
	size_t tEnd = 0;
	char *pLine = new char[1024];
	string strId;
	string strPep;
//	std::regex rgxReverse(m_strRegex);
	string strRev = m_strRegex.substr(1, m_strRegex.npos);
	while(a < tSize)	{
		if(_vs[a].m_vseqBest.empty())	{
			a++;
			continue;
		}
		b = 0;
		while(b < _vs[a].m_vseqBest.size())	{
			tUid = _vs[a].m_vseqBest[b].m_tUid;
			c = 0;
			while(c < _vs[a].m_vseqBest[b].m_vDomains.size())	{
				sprintf(pLine,"%i_%i_%i",(int)_vs[a].m_tId, (int)(b+1), (int)(c+1));
				strId = pLine;
				if(m_mapIdsToPeptides.find(strId) == m_mapIdsToPeptides.end())	{
					cout << strId << "<br />\n";
					c++;
					continue;
				}
				strPep = m_mapIdsToPeptides.find(strId)->second;
				get_pre(_vs[a].m_vseqBest[b].m_strSeq,strPre,_vs[a].m_vseqBest[b].m_vDomains[c].m_lS,_vs[a].m_vseqBest[b].m_vDomains[c].m_lE);
				get_post(_vs[a].m_vseqBest[b].m_strSeq,strPost,_vs[a].m_vseqBest[b].m_vDomains[c].m_lS,_vs[a].m_vseqBest[b].m_vDomains[c].m_lE);
//				if(std::regex_match(_vs[a].m_vseqBest[b].m_strDes,rgxReverse))	{
				if(_vs[a].m_vseqBest[b].m_strDes.find(strRev) == 0)	{
					m_ofOut << "<PeptideEvidence isDecoy=\"true\" post=\"" << strPost;
					m_ofOut << "\" pre=\"" << strPre;
					m_ofOut << "\" end=\"" << (unsigned long)(_vs[a].m_vseqBest[b].m_vDomains[c].m_lE+1);
					m_ofOut << "\" start=\"" << (unsigned long)(_vs[a].m_vseqBest[b].m_vDomains[c].m_lS+1);
					m_ofOut << "\" peptide_ref=\"Pep_" << strPep;
					m_ofOut <<  "\" dBSequence_ref=\"DBSeq" << tUid;
					m_ofOut << "\" id=\"PepEv_" << strId << "\"/>\n";
				}
				else	 {
					m_ofOut << "<PeptideEvidence isDecoy=\"false\" post=\"" << strPost;
					m_ofOut << "\" pre=\"" << strPre;
					m_ofOut << "\" end=\"" << (unsigned long)(_vs[a].m_vseqBest[b].m_vDomains[c].m_lE+1);
					m_ofOut << "\" start=\"" << (unsigned long)(_vs[a].m_vseqBest[b].m_vDomains[c].m_lS+1);
					m_ofOut << "\" peptide_ref=\"Pep_" << strPep;
					m_ofOut <<  "\" dBSequence_ref=\"DBSeq" << tUid;
					m_ofOut << "\" id=\"PepEv_" << strId << "\"/>\n";
				}
				c++;
			}
			b++;
		}
		a++;
	}
	delete pLine;
	return true;
}

// Create the Peptide entries in the mzIndentML file
// Used by the sequence_collection method

bool mzid_report::add_peptides(vector<mspectrum> &_vs)
{

	if(m_ofOut.fail() || !m_ofOut.good())	{
		return false;
	}
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	size_t d = 0;
	size_t tUid = 0;
	string strSeq;
	string strUni;
	string strUniDesc;
	size_t tStart = 0;
	size_t tEnd = 0;
	size_t tS = 0;
	size_t tE = 0;
	map<string,string> mapPeptidesToIds;
	map<string,string>::iterator itPtoI;
	string strOutput;
	string strKey;
	string strId;
	char *pLine = new char[1024];
	while(a < _vs.size())	{
		if(_vs[a].m_vseqBest.empty())	{
			a++;
			continue;
		}
		b = 0;
		while(b < _vs[a].m_vseqBest.size())	{
			tUid = _vs[a].m_vseqBest[b].m_tUid;
			c = 0;
			strOutput = "";
			while(c < _vs[a].m_vseqBest[b].m_vDomains.size())	{
				tStart = _vs[a].m_vseqBest[b].m_vDomains[c].m_lS;
				tEnd = _vs[a].m_vseqBest[b].m_vDomains[c].m_lE;
 				strSeq = _vs[a].m_vseqBest[b].m_strSeq.substr(tStart,tEnd-tStart+1);
				sprintf(pLine,"<Peptide id=\"Pep_");
				strOutput = pLine;
				sprintf(pLine,"%i_%i_%i",(int)_vs[a].m_tId,(int)(b+1),(int)(c+1)); 
				strId = pLine;
				strOutput += pLine;
				sprintf(pLine,"\">\n"); 
				strOutput += pLine;
				strOutput += "\t<PeptideSequence>";
				strKey = strSeq;
				strOutput += strSeq;
				strOutput += "</PeptideSequence>\n";
				d = 0;
				while(d < _vs[a].m_vseqBest[b].m_vDomains[c].m_vAa.size())	{
					strOutput += "\t<Modification monoisotopicMassDelta=\"";
					sprintf(pLine,"%.4lf",_vs[a].m_vseqBest[b].m_vDomains[c].m_vAa[d].m_dMod);
					strKey += " ";
					strKey += pLine;
					strOutput += pLine;
					strOutput += "\" location=\"";
					sprintf(pLine,"%i",_vs[a].m_vseqBest[b].m_vDomains[c].m_vAa[d].m_lPos - _vs[a].m_vseqBest[b].m_vDomains[c].m_lS+1);
					strKey += " ";
					strKey += pLine;
					strOutput += pLine;
					strOutput += "\">\n";
		            get_unimod(_vs[a].m_vseqBest[b].m_vDomains[c].m_vAa[d].m_dMod,strUni,strUniDesc);
					if(!strUni.empty())	{
						strOutput += "\t\t<cvParam accession=\"";
						strOutput += strUni;
						strOutput += "\" cvRef=\"UNIMOD\" name=\"";
						strOutput += strUniDesc;
						strOutput += "\"/>\n";
					}
					strOutput += "\t</Modification>\n";
					d++;
				}
 				strOutput += "</Peptide>\n";

				itPtoI = mapPeptidesToIds.find(strKey);
				if(itPtoI == mapPeptidesToIds.end())	{
					m_ofOut << strOutput;
					m_mapIdsToPeptides[strId] = strId;
					mapPeptidesToIds[strKey] = strId;
				}
				else {
					m_mapIdsToPeptides[strId] = itPtoI->second;
				}
				c++;
			}
			b++;
		}
		a++;
	}
	delete pLine;
	return true;
}

// Parses FASTA description strings into accession numbers and descriptions

bool mzid_report::parse_description(const string &_s,string &_acc,string &_desc)
{
	if(_s.size() == 0)	{
		return false;
	}
	size_t tSpace = _s.find(" ",0);
	if(tSpace == _s.npos)	{
		_acc = _s;
		_desc = "";
		return true;
	}
	_acc = _s.substr(0,tSpace);
	_desc = _s.substr(tSpace+1,_s.length() - tSpace);
	return true;
}

// Obtains the residue N-terminal to the residue indicated by _b

bool mzid_report::get_pre(const string &_s,string &_p,const size_t _b,const size_t _e)
{
	long lStart = (long)_b;
	lStart -= 1;
	_p.erase(_p.begin(),_p.end());
	if(lStart < 0)	{
		lStart = 0;
		_p = '-';
	}
	while(lStart < (long)_b)	{
		_p += _s[lStart];
		lStart++;
	}
	return true;
}

// Obtains the residue C-terminal to the residue indicated by _e

bool mzid_report::get_post(const string &_s,string &_p,const size_t _b,const size_t _e)
{
	size_t tEnd = (long)_e;
	tEnd += 2;
	_p.erase(_p.begin(),_p.end());
	if(tEnd >= _s.size())	{
		tEnd = _s.size();
	}
	size_t a = _e + 1;
	while(a < tEnd)	{
		_p += _s[a];
		a++;
	}
	if(a == _s.size())	{
		_p += '-';
	}
	return true;
}

// Creates a map of modification masses to UniMod identifiers and descriptions
// Called on creation of the mzid_report object: it does not need to be explicitly called

bool mzid_report::load_unimod(void)
{
	int iValue;
	m_mapUnimod.clear();
	m_mapUnimodDescriptions.clear();
	iValue = get_unimod_mass(57.021465);
	m_mapUnimod[iValue] = "UNIMOD:4";
	m_mapUnimodDescriptions[iValue] = "Carbamidomethyl";

	iValue = get_unimod_mass(15.994915);
	m_mapUnimod[iValue] = "UNIMOD:35";
	m_mapUnimodDescriptions[iValue] = "Oxidation";

	iValue = get_unimod_mass(31.989829);
	m_mapUnimod[iValue] = "UNIMOD:425";
	m_mapUnimodDescriptions[iValue] = "Dioxidation";

	iValue = get_unimod_mass(47.984744);
	m_mapUnimod[iValue] = "UNIMOD:435";
	m_mapUnimodDescriptions[iValue] = "Trioxidation";

	iValue = get_unimod_mass(0.9840156);
	m_mapUnimod[iValue] = "UNIMOD:7";
	m_mapUnimodDescriptions[iValue] = "Deamidated";

	iValue = get_unimod_mass(-18.01056);
	m_mapUnimod[iValue] = "UNIMOD:23";
	m_mapUnimodDescriptions[iValue] = "Dehydrated";

	iValue = get_unimod_mass(-17.026548);
	m_mapUnimod[iValue] = "UNIMOD:385";
	m_mapUnimodDescriptions[iValue] = "Ammonia-loss";

	iValue = get_unimod_mass(42.010563);
	m_mapUnimod[iValue] = "UNIMOD:1";
	m_mapUnimodDescriptions[iValue] = "Acetyl";

	iValue = get_unimod_mass(79.966331);
	m_mapUnimod[iValue] = "UNIMOD:21";
	m_mapUnimodDescriptions[iValue] = "Phosphorylation";

	return true;
}

// Returns the UniMod identifier and description string, given a modification mass _m in Daltons

bool mzid_report::get_unimod(double _m,string &_u,string &_d)
{
	int iValue = get_unimod_mass(_m);
	_u.clear();
	_d.clear();
	map<int,string>::iterator itU;
	itU = m_mapUnimod.find(iValue);
	if(itU != m_mapUnimod.end())	{
		_u = itU->second;
		itU = m_mapUnimodDescriptions.find(iValue);
		if(itU != m_mapUnimodDescriptions.end())	{
			_d = itU->second;
		}
	}
	else {
		return false;
	}
	return true;
}

// Returns the UniMod identifier and description string, given a modification mass _m in milliDaltons

bool mzid_report::get_unimod(int _m,string &_u,string &_d)
{
	int iValue = _m;
	_u.clear();
	_d.clear();
	map<int,string>::iterator itU;
	itU = m_mapUnimod.find(iValue);
	if(itU != m_mapUnimod.end())	{
		_u = itU->second;
		itU = m_mapUnimodDescriptions.find(iValue);
		if(itU != m_mapUnimodDescriptions.end())	{
			_d = itU->second;
		}
	}
	else {
		return false;
	}
	return true;
}
