#
# spec file for package cdo
#

Name:           cdo
#BuildRequires:  
Version:        1.5.4
Release:        1
Summary:        Climate Data Operators
License:        GNU GENERAL PUBLIC LICENSE Version 2, June 1991
Group:          Productivity/Graphics/Visualization/Other
Requires:       netcdf
Autoreqprov:    on
URL:            http://code.zmaw.de/projects/cdo
Source0:        cdo-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
CDO is a collection of command line Operators to manipulate and analyse Climate model Data.
Supported data formats are GRIB, netCDF, SERVICE, EXTRA and IEG. There are more than 400
operators available. The following table provides a brief overview of the main categories.


Authors:
--------
    This program was developed at the Max-Planck-Institute for Meteorology.
    Uwe Schulzweida, Uwe.Schulzweida@zmaw.de, is the main author.
    Luis Kornblueh, Luis.Kornblueh@zmaw.de
    Ralf Quast, Ralf.Quast@brockmann-consult.de
    Send questions, comments and bug reports to <http://code.zmaw.de/projects/cdo>


%prep
%setup


%build
./configure --prefix=%{_prefix} --with-netcdf
make 

%install
make DESTDIR=$RPM_BUILD_ROOT install 

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,755)
%doc AUTHORS COPYING ChangeLog NEWS OPERATORS README doc/cdo.pdf doc/cdo_refcard.pdf
%{_prefix}/bin/cdo
#%{_prefix}/bin/cdotest

%changelog -n cdo
* Mon Aug 25 2008 - petri@pik-potsdam.de
- adapted to cdo-1.2.0
* Tue May 20 2008 - petri@pik-potsdam.de
- adapted to cdo-1.1.1
- dont try to include cdotest in the package
* Fri Jan 05 2007 - petri@pik-potsdam.de
- Created initial spec file
