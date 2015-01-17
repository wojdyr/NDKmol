/*  NDKmol - Molecular Viewer on Android NDK
 
 (C) Copyright 2011 - 2012, biochem_fan
 
 This file is part of NDKmol.
 
 NDKmol is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>. */

//  How to compile on Mac OS X
//  g++ -g -std=c++11 -c ../NDKmol/*.cpp
//  g++ -g -std=c++11 -framework GLUT -framework OpenGL -I../NDKmol -o NDKmol *.o NDKmol.cpp
//  ./NDKmol ../3V8T.pdb

// How to compile for WebGL with Emscripten
// ln -s PDB_FILE_YOU_WANT_TO_TEST.pdb initial.pdb
// em++  -std=c++11 -I../NDKmol  -c ../NDKmol/*.cpp
// em++  -std=c++11 -I../NDKmol --preload-file initial.pdb -s TOTAL_MEMORY=100000000 -o NDKmol.html *.o NDKmol.cpp

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <Glut/glut.h>
#else
#define GL_GLEXT_PROTOTYPES
#define EGL_EGLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include "NDKmol/NdkView.h"
#include "NDKmol/Quaternion.h"
#include "NDKmol/Vector3.hpp"

// simple output to stdout for now
#define mesg(fmt, ...) printf(fmt"\n", ##__VA_ARGS__)
#define debug(fmt, ...) printf(fmt"\n", ##__VA_ARGS__)

class Protein;
extern Protein *protein; // global from NDKmol/NdkView.cpp

enum {
  kMenuProteinTrace,
  kMenuProteinThinRibbon,
  kMenuProteinThickRibbon,
  kMenuProteinStrand,
  kMenuProteinTube,
  kMenuProteinBonds,
  kMenuProteinNone,
  kMenuNuclAcidLine,
  kMenuNuclAcidPolygon,
  kMenuLigandSphere,
  kMenuLigandStick,
  kMenuLigandLine,
  kMenuLigandInvisible,
  kMenuColorRainbow,
  kMenuColorChain,
  kMenuColorSS,
  kMenuColorPolarity,
  kMenuColorBfactor,
  kMenuShowMonomer,
  kMenuShowBiological,
  kMenuShowCrystal,
  kMenuToggleSmoothSheets,
  kMenuToggleSideChains,
  kMenuToggleSolvents,
  kMenuToggleUnitCell,
  kMenuToggleHetatmMates,
  kMenuHelp
};

struct WindowState {
  Vector3 obj;
  float cameraZ, slab_near, slab_far;
  Quaternion rotationQ;
  bool fullscreen;
  int normal_width, normal_height; // size of normal (not fullscreen) window

  // values stored when mouse button is pressed
  int mouse_button;
  int start_x, start_y;
  int mouse_modifier; // Shift, Ctrl etc.
  Vector3 current_obj;
  Quaternion currentQ;

  bool menu_in_use;
  int protein_mode;
  int nucleic_acid_mode;
  int hetatm_mode;
  int symmetry_mode;
  int color_mode;
  bool show_sidechains;
  bool show_unitcell;
  bool show_solvents;
  bool smoothen;
  bool symop_hetatms;
};

static WindowState w;

static void init_state() {
  w.rotationQ.w = -1;
  w.fullscreen = false;
  w.normal_width = 800;
  w.normal_height = 600;
  w.mouse_button = -1;
  w.menu_in_use = false;
  w.protein_mode = MAINCHAIN_THICKRIBBON;
  w.nucleic_acid_mode = BASE_LINE;
  w.hetatm_mode = HETATM_STICK;
  w.symmetry_mode = SYMOP_BIOMT;
  w.color_mode = COLOR_CHAINBOW;
  w.show_sidechains = false;
  w.show_unitcell = false;
  w.show_solvents = false;
  w.smoothen = true;
  w.symop_hetatms = false;
}

static Vector3 operator+(const Vector3& a, const Vector3& b) {
  return Vector3(a.x+b.x, a.y+b.y, a.z+b.z);
}

static void render() {
  glClear(GL_COLOR_BUFFER_BIT);
  float ax, ay, az;
  w.rotationQ.getAxis(&ax, &ay, &az);
  nativeSetScene(w.obj.x, w.obj.y, w.obj.z, ax, ay, az,
                 w.rotationQ.getAngle(), w.cameraZ, w.slab_near, w.slab_far);
  nativeGLRender();
  glutSwapBuffers();
}

static void rebuild_scene() {
  bool reset_view = false;
  buildScene(w.protein_mode, w.hetatm_mode, w.symmetry_mode, w.color_mode,
             w.show_sidechains, w.show_unitcell,
             w.nucleic_acid_mode, w.show_solvents,
             reset_view, !w.smoothen, w.symop_hetatms);
  glutPostRedisplay();
}

static void on_change_size(int w, int h) {
  glViewport(0, 0, w, h);
  nativeGLResize(w, h);
}

static void toggle_fullscreen() {
  w.fullscreen = !w.fullscreen;
  if (w.fullscreen) {
    w.normal_width = glutGet(GLUT_WINDOW_WIDTH);
    w.normal_height = glutGet(GLUT_WINDOW_HEIGHT);
    glutFullScreen();
  } else {
    glutReshapeWindow(w.normal_width, w.normal_height);
  }
}

static void show_help() {
  mesg("Help is not implemented yet.");
}

static void menu_handler(int option) {
  switch (option) {
    case kMenuProteinTrace:
      w.protein_mode = MAINCHAIN_TRACE;
      break;
    case kMenuProteinThinRibbon:
      w.protein_mode = MAINCHAIN_THINRIBBON;
      break;
    case kMenuProteinThickRibbon:
      w.protein_mode = MAINCHAIN_THICKRIBBON;
      break;
    case kMenuProteinStrand:
      w.protein_mode = MAINCHAIN_STRAND;
      break;
    case kMenuProteinTube:
      w.protein_mode = MAINCHAIN_TUBE;
      break;
    case kMenuProteinBonds:
      w.protein_mode = MAINCHAIN_BONDS;
      break;
    case kMenuProteinNone:
      w.protein_mode = MAINCHAIN_NONE;
      break;
    case kMenuNuclAcidLine:
      w.nucleic_acid_mode = (w.nucleic_acid_mode == BASE_LINE ? BASE_NONE
                                                              : BASE_LINE);
      break;
    case kMenuNuclAcidPolygon:
      w.nucleic_acid_mode = (w.nucleic_acid_mode == BASE_POLYGON ? BASE_NONE
                                                              : BASE_POLYGON);
      break;
    case kMenuLigandSphere:
      w.hetatm_mode = HETATM_SPHERE;
      break;
    case kMenuLigandStick:
      w.hetatm_mode = HETATM_STICK;
      break;
    case kMenuLigandLine:
      w.hetatm_mode = HETATM_LINE;
      break;
    case kMenuLigandInvisible:
      w.hetatm_mode = HETATM_NONE;
      break;
    case kMenuColorRainbow:
      w.color_mode = COLOR_CHAINBOW;
      break;
    case kMenuColorChain:
      w.color_mode = COLOR_CHAIN;
      break;
    case kMenuColorSS:
      w.color_mode = COLOR_SS;
      break;
    case kMenuColorPolarity:
      w.color_mode = COLOR_POLAR;
      break;
    case kMenuColorBfactor:
      w.color_mode = COLOR_B_FACTOR;
      break;
    case kMenuShowMonomer:
      w.symmetry_mode = SYMOP_NONE;
      break;
    case kMenuShowBiological:
      w.symmetry_mode = SYMOP_BIOMT;
      break;
    case kMenuShowCrystal:
      w.symmetry_mode = SYMOP_PACKING;
      break;
    case kMenuToggleSmoothSheets:
      w.smoothen = !w.smoothen;
      break;
    case kMenuToggleSideChains:
      w.show_sidechains = !w.show_sidechains;
      break;
    case kMenuToggleUnitCell:
      w.show_unitcell = !w.show_unitcell;
      break;
    case kMenuToggleSolvents:
      w.show_solvents = !w.show_solvents;
      break;
    case kMenuToggleHetatmMates:
      w.symop_hetatms = !w.symop_hetatms;
      break;
    case kMenuHelp:
      show_help();
      return; // skip rebuild_scene()
  }
  rebuild_scene();
}

static Vector3 pan_by(float dx, float dy) {
  float f = -w.cameraZ * 0.0005;
  Vector3 v(f * dx, -f * dy, 0);
  w.rotationQ.rotateVector(v.x, v.y, v.z, &v.x, &v.y, &v.z);
  return v;
}

static void on_key(unsigned char key, int /*x*/, int /*y*/) {
  switch(key) {
    case '+':
    case '=':
      w.cameraZ /= 1.2;
      glutPostRedisplay();
      break;
    case '-':
    case '_':
      w.cameraZ *= 1.2;
      glutPostRedisplay();
      break;
    case 'f':
      toggle_fullscreen();
      break;
    case 's': menu_handler(kMenuToggleSmoothSheets); break;
    case 'd': menu_handler(kMenuToggleSideChains); break;
    case 'c': menu_handler(kMenuToggleUnitCell); break;
    case 'w': menu_handler(kMenuToggleSolvents); break;
    case 'm': menu_handler(kMenuToggleHetatmMates); break;
    case '1': menu_handler(kMenuColorRainbow); break;
    case '2': menu_handler(kMenuColorChain); break;
    case '3': menu_handler(kMenuColorSS); break;
    case '4': menu_handler(kMenuColorPolarity); break;
    case '5': menu_handler(kMenuColorBfactor); break;
    case '7': menu_handler(kMenuShowMonomer); break;
    case '8': menu_handler(kMenuShowBiological); break;
    case '9': menu_handler(kMenuShowCrystal); break;
    case 'b': menu_handler(kMenuNuclAcidPolygon); break;
    case 'n': menu_handler(kMenuNuclAcidLine); break;
    case 'y': menu_handler(kMenuLigandSphere); break;
    case 'u': menu_handler(kMenuLigandStick); break;
    case 'i': menu_handler(kMenuLigandLine); break;
    case 'o': menu_handler(kMenuLigandInvisible); break;

    case 27: // ESC
    case 'q':
      exit(0);
      break;
    default:
      mesg("key '%c' does nothing.", key);
  }
}

static void on_special_key(int key, int /*x*/, int /*y*/) {
  switch (key) {
    case GLUT_KEY_F1: menu_handler(kMenuProteinTrace); break;
    case GLUT_KEY_F2: menu_handler(kMenuProteinThinRibbon); break;
    case GLUT_KEY_F3: menu_handler(kMenuProteinThickRibbon); break;
    case GLUT_KEY_F4: menu_handler(kMenuProteinStrand); break;
    case GLUT_KEY_F5: menu_handler(kMenuProteinTube); break;
    case GLUT_KEY_F6: menu_handler(kMenuProteinBonds); break;
    case GLUT_KEY_F7: menu_handler(kMenuProteinNone); break;
    case GLUT_KEY_F8: break;
    case GLUT_KEY_F9: break;
    case GLUT_KEY_F10: break;
    case GLUT_KEY_F11: toggle_fullscreen(); break;
    case GLUT_KEY_F12: break;

    case GLUT_KEY_LEFT:
      w.obj = w.obj + pan_by(-40, 0);
      rebuild_scene();
      break;
    case GLUT_KEY_RIGHT:
      w.obj = w.obj + pan_by(+40, 0);
      rebuild_scene();
      break;
    case GLUT_KEY_UP:
      w.obj = w.obj + pan_by(0, -40);
      rebuild_scene();
      break;
    case GLUT_KEY_DOWN:
      w.obj = w.obj + pan_by(0, +40);
      rebuild_scene();
      break;
    case GLUT_KEY_PAGE_UP: break;
    case GLUT_KEY_PAGE_DOWN: break;
    case GLUT_KEY_HOME: break;
    case GLUT_KEY_END: break;
    case GLUT_KEY_INSERT: break;
    default:
      break;
  }
}

static void on_mouse_move(int x, int y) {
  if (w.mouse_button == GLUT_LEFT_BUTTON) {
    float dx = (x - w.start_x) / (float)glutGet(GLUT_WINDOW_WIDTH);
    float dy = (y - w.start_y) / (float)glutGet(GLUT_WINDOW_HEIGHT);
    float r = sqrt(dx * dx + dy * dy);
    if (r == 0)
        return;
    float rs = sin(r * M_PI) / r;
    Quaternion dq(rs * dy, rs * dx, 0, cos(r * M_PI));
    w.rotationQ = Quaternion::multiply(dq, w.currentQ);
    glutPostRedisplay();
  }
  else if (w.mouse_button == GLUT_MIDDLE_BUTTON) {
    w.obj = w.current_obj + pan_by(x - w.start_x, y - w.start_y);
    glutPostRedisplay();
  }
}

static void on_mouse_button(int button, int state, int x, int y) {
  w.mouse_modifier = glutGetModifiers();
  if (button == GLUT_LEFT_BUTTON || button == GLUT_MIDDLE_BUTTON) {
    if (state == GLUT_DOWN) {
      w.mouse_button = button;
      w.start_x = x;
      w.start_y = y;
      w.current_obj = w.obj;
      w.currentQ = w.rotationQ;
      // TODO: with Shift+Left: zooming + rotation like in JMol
    } else {
      on_mouse_move(x, y);
      w.mouse_button = -1;
    }
  } else if ((button == 3 || button == 4)&& state == GLUT_DOWN) { // scroll
    float zoom_factor = 1.1f;
    if (w.mouse_modifier & GLUT_ACTIVE_SHIFT)
      zoom_factor = 1.01f;
    if (button == 3)
      w.cameraZ /= zoom_factor;
    else // button == 4
      w.cameraZ *= zoom_factor;
    glutPostRedisplay();
  }
}

static void on_menu_status(int status, int /*x*/, int /*y*/) {
  w.menu_in_use = (status == GLUT_MENU_IN_USE);
}

static void create_menu() {
  int polymer_menu = glutCreateMenu(menu_handler);
  glutAddMenuEntry("mainchain trace [F1]", kMenuProteinTrace);
  glutAddMenuEntry("thin ribbon [F2]", kMenuProteinThinRibbon);
  glutAddMenuEntry("thick ribbon [F3]", kMenuProteinThickRibbon);
  glutAddMenuEntry("strand [F4]", kMenuProteinStrand);
  glutAddMenuEntry("B-factor tube [F5]", kMenuProteinTube);
  glutAddMenuEntry("all bonds [F6]", kMenuProteinBonds);
  glutAddMenuEntry("invisible [F7]", kMenuProteinNone);

  int ligands_menu = glutCreateMenu(menu_handler);
  glutAddMenuEntry("Sphere [y]", kMenuLigandSphere);
  glutAddMenuEntry("Stick [u]", kMenuLigandStick);
  glutAddMenuEntry("Line [i]", kMenuLigandLine);
  glutAddMenuEntry("Invisible [o]", kMenuLigandInvisible);

  int color_menu = glutCreateMenu(menu_handler);
  glutAddMenuEntry("Rainbow [1]", kMenuColorRainbow);
  glutAddMenuEntry("Chain [2]", kMenuColorChain);
  glutAddMenuEntry("Structure [3]", kMenuColorSS);
  glutAddMenuEntry("Polarity [4]", kMenuColorPolarity);
  glutAddMenuEntry("B-factor [5]", kMenuColorBfactor);

  int packing_menu = glutCreateMenu(menu_handler);
  glutAddMenuEntry("Monomer [7]", kMenuShowMonomer);
  glutAddMenuEntry("Biological Unit [8]", kMenuShowBiological);
  glutAddMenuEntry("Crystal Packing [9]", kMenuShowCrystal);

  int toggle_menu = glutCreateMenu(menu_handler);
  glutAddMenuEntry("Smooth Sheets [s]", kMenuToggleSmoothSheets);
  glutAddMenuEntry("Sidechains [d]", kMenuToggleSideChains);
  glutAddMenuEntry("Solvents [w]", kMenuToggleSolvents);
  glutAddMenuEntry("Unit Cell [c]", kMenuToggleUnitCell);
  glutAddMenuEntry("Nucleic Acid Polygons [b]", kMenuNuclAcidPolygon);
  glutAddMenuEntry("Nucleic Acid Lines [n]", kMenuNuclAcidLine);
  glutAddMenuEntry("HETATM symmetry mates [m]", kMenuToggleHetatmMates);


  glutCreateMenu(menu_handler); // main menu
  glutAddSubMenu("Polymer as", polymer_menu);
  glutAddSubMenu("Ligands as", ligands_menu);
  glutAddSubMenu("Color by", color_menu);
  glutAddSubMenu("Show", packing_menu);
  glutAddSubMenu("Toggle", toggle_menu);
  glutAddMenuEntry("Help", kMenuHelp);

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// returns uppercased extension
static std::string get_extension(const char* filename) {
  std::string ext;
  const char* last_dot = strrchr(filename, '.');
  if (last_dot && last_dot > filename && strcmp(last_dot-1, ".gz") == 0)
    last_dot = strrchr(last_dot-1, '.');
  if (last_dot)
    for (int i = 0; last_dot[i+1] != '\0'; ++i)
      ext.append(1, toupper(last_dot[i+1]));
  return ext;
}

static void open_file(const char* filename) {
  std::string ext = get_extension(filename);
  if (ext == "PDB" || ext == "ENT")
    nativeLoadProtein(filename);
  else if (ext == "SDF" || ext == "MOL")
    nativeLoadSDF(filename);
  else if (ext == "CCP4" || ext == "CCP4.GZ" || ext == "MAP")
    nativeLoadCCP4(filename);
  else
    mesg("Unrecognized filetype (%s) of: %s", ext.c_str(), filename);
}

int main(int argc, char **argv) {
  init_state();

  const char *filename = "res/raw/initial.pdb"; //FIXME
  std::string title = "NDKmol";
  if (argc > 1) {
    filename = argv[1];
    title = std::string("NDKmol - ") + filename;
  }

  glutInit(&argc, argv);
  glutInitWindowSize(w.normal_width, w.normal_height);
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow(title.c_str());
  glutDisplayFunc(render);
  glutReshapeFunc(on_change_size);
  glutKeyboardFunc(on_key);
  glutSpecialFunc(on_special_key);
  glutMouseFunc(on_mouse_button);
  glutMotionFunc(on_mouse_move);
  glutMenuStatusFunc(on_menu_status);

  open_file(filename); //TODO: handling errors
  if (protein != NULL && argc > 2 &&
      get_extension(argv[2]).substr(0,4) == "CCP4")
    nativeLoadCCP4(argv[2]);

  nativeAdjustZoom(&w.obj.x, &w.obj.y, &w.obj.z,
                   &w.cameraZ, &w.slab_near, &w.slab_far, false);
  rebuild_scene();
  nativeGLInit();

  create_menu();
  glutMainLoop();

  return 0;
}

// vim: et:ts=2:sw=2
