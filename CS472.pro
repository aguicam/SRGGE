TEMPLATE = app
TARGET = ViewerPro
QT += widgets opengl 
CONFIG += qt debug_and_release
RESOURCES   = CS472.qrc

Release:OBJECTS_DIR = release/.obj
Release:MOC_DIR     = release/.moc
Debug:OBJECTS_DIR   = debug/.obj
Debug:MOC_DIR       = debug/.moc

INCLUDEPATH     += . geometry camera lighting Eigen GL glm

win32-msvc2015{
	RC_FILE        += CS472.rc
        LIBS           += -lopengl32
	QMAKE_CXXFLAGS += /MP /Zi
	QMAKE_LFLAGS   += /MACHINE:X64
}

macx{
        QMAKE_SONAME_PREFIX = @executable_path/../Frameworks
        QMAKE_LFLAGS += -Wl,-rpath,@executable_path/../Frameworks
	QMAKE_LFLAGS   -= -mmacosx-version-min=10.8
	QMAKE_LFLAGS   += -mmacosx-version-min=10.9
	QMAKE_CXXFLAGS -= -mmacosx-version-min=10.8
	QMAKE_CXXFLAGS += -mmacosx-version-min=10.9
	QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9
	ICON = CS472.icns
}

unix:!macx{
        CONFIG += C++11
}

mingw{
        LIBS += -lopengl32
}

# Input
HEADERS +=	MainWindow.h		\
		HW.h			\
		dummy.h			\
		camera/Camera.h \
                Helpers/mesh_io.h \
                Helpers/triangle_mesh.h \
    Helpers/mesh_importer.h \
    Helpers/vis.h \
    SRGGE/Viewer1.h \
    SRGGE/Viewer2.h \
    SRGGE/viewer3.h \
    SRGGE/viewer4.h \
    SRGGE/Viewer5.h \
    Helpers/stb_image.h
		
SOURCES +=	main.cpp		\ 
		MainWindow.cpp 		\
		HW.cpp			\
		dummy.cpp		\
		camera/Camera.cpp \
                Helpers/mesh_io.cc \
                Helpers/triangle_mesh.cc \
    Helpers/mesh_importer.cpp \
    Helpers/vis.cpp \
    SRGGE/Viewer1.cpp \
    SRGGE/Viewer2.cpp \
    SRGGE/Viewer3.cpp \
    SRGGE/Viewer4.cpp \
    SRGGE/Viewer5.cpp

