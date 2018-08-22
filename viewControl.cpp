#include "viewControl.h"
#include <QPainter>
#include <QPrintDialog>
#include <algorithm>
#include <array>
#include <iostream>

namespace Ui {
    class MainWindow;
}

inline bool disPoint(QPoint a, QPoint b, int c){
    return (a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y()) <=c*c;
}

inline float distPoint(QPoint a, QPoint b)
{
    return (a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y());
}

ViewControl::ViewControl(QFrame* frame, Viewer* viewer, MainWindow* window): QWidget(frame) {
    parentViewer = viewer;
    parentWindow = window;
    image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
    image.fill(qRgba(255, 255, 255, 1));

    penColor = Qt::black;
    penWidth = 4;
    skeWidth = 4;

    setWindowFlags(windowFlags() | Qt::FramelessWindowHint);
    setAttribute(Qt::WA_TranslucentBackground);

    modified = false;
    isDrawing = false;

    leftmax = upmax = VIEW_SIZE;
    rightmax = downmax = 0;
    skerightmax = skedownmax = 0;
    skeleftmax = skeupmax = VIEW_SIZE;
}

void ViewControl::changeImage(QImage newImage){
    image = newImage;
}

void ViewControl::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    if(isDrawing) {
        painter.drawImage(0,0,tempImage);
    } else {
        painter.drawImage(0,0,image);
    }
}

void ViewControl::keyPressEvent(QKeyEvent *event){
    if (event->key()==Qt::Key_1){
        //TODO
    } else if (event->key()==Qt::Key_9){
    //TODO
        if (parentWindow->inner_mode==SKETC_MODE){
            // undo
            int nsketches = sketchesbuffer.size();
            if (nsketches>0&&cursketch>=0&&cursketch<nsketches) {
                sketches = sketchesbuffer[cursketch];
                image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
                image.fill(qRgba(255, 255, 255, 1));
                paint(image);
            }
        }
    }
}

void ViewControl::mousePressEvent(QMouseEvent *event) {

    parentWindow->moveViewControl();

    if(parentWindow->inner_mode == REFINE_SELECT_MODE || parentWindow->inner_mode == REFINE_DEFORM_MODE){
        drawPoints.erase(drawPoints.begin(), drawPoints.end());
        parentViewer->removeAfterCurrentpos();
        if(parentWindow->inner_mode == REFINE_SELECT_MODE) {
            if (event->button() == Qt::LeftButton){
                penWidth = 10;
                penColor = Qt::red;
                //projectSelection(true);
                //selectHandle();
            } else if (event->button() == Qt::RightButton){
                penWidth = 20;
                penColor = Qt::gray;
                //projectSelection(false);
               // eraseHandle();
            }
            lineStartPoint = startPoint = event->pos();
            drawPoints.push_back(QPoint(event->pos().x(), VIEW_SIZE - event->pos().y()));
            isDrawing = true;
            setMaxBound(startPoint);

        } else if (event->button() == Qt::LeftButton && parentWindow->inner_mode == REFINE_DEFORM_MODE){
            penColor = Qt::blue;
            lineStartPoint = startPoint = event->pos();
            isDrawing = true;
            setMaxBound(startPoint);
        }
    }
    else if(parentWindow->inner_mode == LEFT_LABEL || parentWindow->inner_mode == RIGHT_LABEL)
    {
        //TODO
         if(event->button() == Qt::LeftButton)
         {
            isdrawdeformline = true;
            sidedeformline.clear();
            sidedeformline.push_back(event->pos());
            isdrawline = true;

         }
    }
    else if (parentWindow->inner_mode == SKETCH_FINE_MODE){
        //TODO
         if(event->button() == Qt::LeftButton)
         {
             if(is_frontmode)
             {
                 if(erasesketchid<0)
                 {
                     isdrawwrinkle = true;
                     isdrawdeformline = false;
                     wrinkleline.clear();
                     wrinkleline.push_back(event->pos());
                 }
                 else
                 {
                     isdrawdeformline = true;
                     isdrawwrinkle = false;
                     deformline.clear();
                     deformline.push_back(event->pos());
                 }
                 isdrawline = true;
             }
             else
             {
                /* if(sideerasesketchid<0)
                 {
                    // isdrawwrinkle = true;
                     isdrawdeformline = false;
                     sidewrinkleline.clear();
                     sidewrinkleline.push_back(event->pos());
                 }
                 else*/
                 {
                     isdrawdeformline = true;
                    // isdrawwrinkle = false;
                     sidedeformline.clear();
                     //sidedeformline.push_back(event->pos());
                     sidedeformline.push_back(QPoint(event->pos().x(), VIEW_SIZE-event->pos().y()));
                 }
                 isdrawline = true;
             }
         }

         if(event->button() == Qt::RightButton)
         {
             if(is_frontmode)
             {
                 iseraseline = true;
                 eraseline.clear();
                 eraseline.push_back(event->pos());
                 erasedPoints.clear();
                 erasedColor.clear();
             }
             else
             {
                /* iseraseline = true;
                 sideeraseline.clear();
                 sideeraseline.push_back(event->pos());
                 sideerasedColor.clear();
                 sideerasedPoints.clear();*/
             }

         }
    }

}

void ViewControl::mouseMoveEvent(QMouseEvent *event) {
    if(parentWindow->inner_mode == REFINE_SELECT_MODE || parentWindow->inner_mode == REFINE_DEFORM_MODE){
        if(parentWindow->inner_mode == REFINE_SELECT_MODE) {
            if (event->buttons()&Qt::LeftButton){
                endPoint = event->pos();
                setMaxBound(endPoint);
                drawPoints.push_back(QPoint(event->pos().x(), VIEW_SIZE - event->pos().y()));
                tempImage = image;
                paint(image);
                //projectSelection(true);
                selectHandle(event->pos());

            } else if (event->buttons()&Qt::RightButton){
                endPoint = event->pos();
                setMaxBound(endPoint);
                drawPoints.push_back(QPoint(event->pos().x(), VIEW_SIZE - event->pos().y()));
                tempImage = image;
                paint(image);
                //projectSelection(false);
                eraseHandle(event->pos());

            }
        } else if((event->buttons()&Qt::LeftButton) && parentWindow->inner_mode == REFINE_DEFORM_MODE) {
            drawPoints.push_back(QPoint(event->pos().x(), VIEW_SIZE - event->pos().y()));
            setMaxBound(endPoint);
            endPoint = event->pos();
            tempImage = image;
            paint(image);
        }
    }
    else if(parentWindow->inner_mode == LEFT_LABEL || parentWindow->inner_mode == RIGHT_LABEL)
    {

        if(isdrawdeformline)
        {
            sidedeformline.push_back(event->pos());

            tempImage = image;
            paint(image);
        }
    }
    else if (parentWindow->inner_mode == SKETCH_FINE_MODE){
        //TODO
        if((event->buttons()&Qt::LeftButton) && isdrawline) {
            if(is_frontmode)
            {
                if(isdrawdeformline)
                {
                    deformline.push_back(event->pos());
                }
                else if(isdrawwrinkle)
                {
                    wrinkleline.push_back(event->pos());
                }

            }
            else
            {
                if(isdrawdeformline)
                {
                    sidedeformline.push_back(event->pos());
                }
              /*  else if(isdrawwrinkle)
                {
                    sidewrinkleline.push_back(event->pos());
                }*/
            }

             tempImage = image;
             paint(image);
        }

        if((event->buttons()&Qt::RightButton) && iseraseline)
        {
            if(is_frontmode)
            {
                eraseline.push_back(event->pos());
                temporalEraseSketch(event->pos());
            }
            else
            {
                sideeraseline.push_back(event->pos());
                temporalEraseSideSketch(event->pos());
            }

           tempImage = image;
           paint(image);
        }
    }
}

void ViewControl::mouseReleaseEvent(QMouseEvent *event) {
    if(parentWindow->inner_mode == REFINE_SELECT_MODE || parentWindow->inner_mode == REFINE_DEFORM_MODE){
        endPoint = event->pos();
        setMaxBound(endPoint);
        isDrawing = false;
        finalImage = image;
        drawPoints.push_back(QPoint(event->pos().x(), VIEW_SIZE - event->pos().y()));
        if(event->button() == Qt::LeftButton && parentWindow->inner_mode == REFINE_SELECT_MODE){
           // projectSelection(true);
            selectHandle(event->pos());
        } else if (event->button() == Qt::RightButton && parentWindow->inner_mode == REFINE_SELECT_MODE){
            //projectSelection(false);
            eraseHandle(event->pos());
        } else if (event->button() == Qt::LeftButton && parentWindow->inner_mode == REFINE_DEFORM_MODE){
            if(drawPoints.size() > 5) {
                sample_cnt = 0;
                saveCroppedImage();
                int gesture_num = parentWindow->run_gesture();
                if (gesture_num <= 10){
                    parentViewer->sampled = false;
                    printf("drag num %d\n", gesture_num);
                   // doShiftArea((gesture_num<=5), gesture_num%5?gesture_num%5:5);
                    doLocalDeform((gesture_num<=5), gesture_num%5?gesture_num%5:5);
                    shifted = true;
                } else if (gesture_num == 11||gesture_num == 12){                                
                    parentViewer->removeAfterCurrentpos();
                    QPoint a = lineStartPoint;
                    QPoint b = endPoint;
                    cout<<a.x()<<" "<<a.y()<<" "<<b.x()<<" "<<b.y()<<endl;
                    if(distPoint(a, b)<400)
                    {
                        printf("selectArea\n");
                       // if (shifted) doClear();
                       // parentViewer->removeAfterCurrentpos();
                        projectNeighborSelection();
                    }
                    else
                    {
                        if(parentViewer->isROISelection)
                        {
                            doLineDeform(a, b);
                        } else if(!parentViewer->isHandleSelection){
                            printf("silhouette Deform\n");
                            // TODO
                        }

                    }

                    //parentViewer->upsampleArea(leftmax, rightmax, upmax, downmax, drawPoints);
    //                selectedLine.clear();
    //                selectedLine = drawPoints;

    //                parentViewer->upsamplelinearea(leftmax, rightmax, upmax, downmax, drawPoints);
                    //selectLine();
                    //doShiftLinearea(1,1);
    //                lineSelected = true;
                }else if (gesture_num == 13){
                    if (parentViewer->sampled){
                        parentWindow->undo();
                        parentViewer->sampled = false;
                    }
                    printf("UNDO\n");
                    parentWindow->undo();
                } else if (gesture_num == 14){
                    if (parentViewer->sampled){
                        parentWindow->redo();
                        parentViewer->sampled = false;
                    }
                    printf("REDO\n");
                    parentWindow->redo();
                } else {
                    printf("Meaningless Gesture.\n");
                }
            }

        }
        penWidth = 4;

        image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
        image.fill(qRgba(255, 255, 255, 1));
        finalImage = image;
        leftmax = upmax = VIEW_SIZE;
        rightmax = downmax = 0;
        endPoint = startPoint = QPoint(-1, -1);
        update();
    }
    else if(parentWindow->inner_mode == LEFT_LABEL || parentWindow->inner_mode == RIGHT_LABEL)
    {
        image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
        image.fill(qRgba(255, 255, 255, 1));
        if(isdrawdeformline)
        {
            sidedeformline.push_back(event->pos());
            if(sidedeformline.size()>3)
            {
                mapSideDeformline();
            }
            sidedeformline.clear();
            isdrawdeformline = false;

            parentWindow->run_sidedeform();
        }

        paint(image);
        isdrawline = false;
    }
    else if (parentWindow->inner_mode == SKETCH_FINE_MODE){
        //TODO
        if(event->button() == Qt::LeftButton && isdrawline){
            image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
            image.fill(qRgba(255, 255, 255, 1));
            if(is_frontmode)
            {
                if(isdrawwrinkle)
                {
                    wrinkleline.push_back(event->pos());
                    wrinkles.push_back(wrinkleline);
                    wrinkleline.clear();
                    isdrawwrinkle = false;
                }
                else if(isdrawdeformline)
                {
                    deformline.push_back(event->pos());
                    resampleline(deformline);
                    if(deformline.size()>3)
                    {
                        mappingDeformline(erasesketchid);
                    }
                    deformline.clear();
                    isdrawdeformline = false;
                }
            }
            else
            {
               /* if(isdrawwrinkle)
                {
                    sidewrinkleline.push_back(event->pos());
                    sidewrinkles.push_back(sidewrinkleline);
                    sidewrinkleline.clear();
                    isdrawwrinkle = false;
                }
                else if(isdrawdeformline)*/
                if(isdrawdeformline)
                {
                    sidedeformline.push_back(event->pos());
                   // resampleline(sidedeformline);
                    if(sidedeformline.size()>3)
                    {
                        mapSideDeformline();
                    }
                    sidedeformline.clear();
                    isdrawdeformline = false;
                }
               // is_sideedit = true;
            }

            paint(image);
            paintfrontimg();
            parentWindow->finePanel->setImage(frontimg);
            isdrawline = false;
        }

        if(event->button() == Qt::RightButton && iseraseline)
        {
           image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
           image.fill(qRgba(255, 255, 255, 1));

           if(is_frontmode)
           {
               eraseline.push_back(event->pos());

               for(int i=0; i<erasedPoints.size(); i++)
               {
                   int ii = erasedPoints[i].first;
                   int jj = erasedPoints[i].second;
                   sketchescolor[ii][jj] = skeColor;
               }
               erasedPoints.clear(); erasedColor.clear();
               eraseSketch();
               eraseline.clear();
           }
          /* else
           {
               sideeraseline.push_back(event->pos());

               for(int i=0; i<sideerasedPoints.size(); i++)
               {
                   int ii = sideerasedPoints[i].first;
                   int jj = sideerasedPoints[i].second;
                   sidesketchescolor[ii][jj] = skeColor;
               }
               sideerasedPoints.clear(); sideerasedColor.clear();
               eraseSideSketch();
               sideeraseline.clear();
               is_sideedit = true;
           }*/

           paint(image);
           iseraseline = false;
        }
        endPoint = startPoint = QPoint(-1, -1);
    }

}

void ViewControl::clearSelection()
{
   int nq = parentViewer->hightlight_deformneighbor_rawindex.size();
   for(int i=0; i<nq; i++)
   {
       parentViewer->hightlight_deformneighbor_rawindex[i] = false;
       parentViewer->hightlight_selection_rawindex[i] = false;
   }
   parentViewer->isROISelection = false;
   parentViewer->isHandleSelection = false;
   parentViewer->updateGL();
}

void ViewControl::setMaxBound(const QPoint & point){
    if (point.x() < 0 || point.y() < 0) return;
    if (point.x() < leftmax) leftmax = point.x();
    if (point.x() > rightmax) rightmax = point.x();
    if (point.y() < upmax) upmax = point.y();
    if (point.y() > downmax) downmax = point.y();
}

void ViewControl::setSideSkeMaxBound(const QPoint & point)
{
    if (point.x() < 0 || point.y() < 0) return;
    if (point.x() < sideskeleftmax) sideskeleftmax = point.x();
    if (point.x() > sideskerightmax) sideskerightmax = point.x();
    if (point.y() < sideskeupmax) sideskeupmax = point.y();
    if (point.y() > sideskedownmax) sideskedownmax = point.y();
}

void ViewControl::setSkeMaxBound(const QPoint & point){
    if (point.x() < 0 || point.y() < 0) return;
    if (point.x() < skeleftmax) skeleftmax = point.x();
    if (point.x() > skerightmax) skerightmax = point.x();
    if (point.y() < skeupmax) skeupmax = point.y();
    if (point.y() > skedownmax) skedownmax = point.y();
}

void ViewControl::saveCroppedImage(){
    int boundsize = std::max(rightmax-leftmax, downmax-upmax);
    boundsize *= 1.4;
    int centerpointx = (rightmax+leftmax)/2;
    int centerpointy = (downmax+upmax)/2;
    int newleft = centerpointx - boundsize/2;
    int newup = centerpointy - boundsize/2;
    QImage visibleImage = finalImage.copy(newleft, newup, boundsize, boundsize).convertToFormat(QImage::Format_RGB32, Qt::ColorOnly);
    int iw = visibleImage.size().width(), ih = visibleImage.size().height();
    for (int x = 0; x < iw; ++x) {
        for (int y = 0; y < ih; ++y) {
          if (visibleImage.pixelColor(x, y).red() > 252) continue; //white
          else if (visibleImage.pixelColor(x, y).blue() < 250) visibleImage.setPixel(x, y, qRgb(255, 255, 255)); // black -> white
          else visibleImage.setPixel(x, y, qRgb(0, 0, 0)); // red -> black
        }
    }
    visibleImage.save(std::string(localfolder()+"\\models\\front\\gesture.png").c_str(), "png");
}

void ViewControl::saveSketchImage(){
  frontimg.save(std::string(localfolder()+"\\models\\front\\saved.png").c_str(), "png");
}

void ViewControl::saveUncroppedImage(){
    //pass
}

void ViewControl::projectAllPoints(){
    int currentpos = parentViewer->currentpos;
    int nq = parentViewer->raw_vertices[currentpos].size();
    rawproject.clear();

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());


    //const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        //float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        //if(cur_angle < 0){ // the front points
            gluProject(parentViewer->raw_vertices[currentpos][i].x(), parentViewer->raw_vertices[currentpos][i].y(), parentViewer->raw_vertices[currentpos][i].z(),
                modelview.data(), projection.data(), viewport.data(),
                &current_x, &current_y, &current_z);
            rawproject.push_back(make_pair(current_x, VIEW_SIZE - current_y));
        //}
    }
}

void ViewControl::selectHandle(QPoint& p)
{
    int currentpos = parentViewer->currentpos;

    // start finding
    int nq = rawproject.size();

    area_x_left = p.x()-20;  area_x_right = p.x()+20;
    area_y_up = p.y()+20;  area_y_down = p.y()-20;

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            int int_current_x = rawproject[i].first, int_current_y = rawproject[i].second;
            QPoint a (int_current_x, int_current_y);
            if(disPoint(a, p, 5))
                 parentViewer->hightlight_selection_rawindex[i] = true;
        }
    }

    parentViewer->isHandleSelection = true;
    parentViewer->updateGL();

}

void ViewControl::eraseHandle(QPoint& p)
{
    int currentpos = parentViewer->currentpos;

    // start finding
    int nq = rawproject.size();

    area_x_left = p.x()-20;  area_x_right = p.x()+20;
    area_y_up = p.y()+20;  area_y_down = p.y()-20;

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            int int_current_x = rawproject[i].first, int_current_y = rawproject[i].second;
            QPoint a (int_current_x, int_current_y);
            if(disPoint(a, p, 10))
                 parentViewer->hightlight_selection_rawindex[i] = false;
        }
    }

    parentViewer->updateGL();
}

void ViewControl::projectSelection(bool draw){
    int currentpos = parentViewer->currentpos;
    shift_candidates.clear();

    QImage shiftImage(std::move(finalImage.convertToFormat(QImage::Format_RGB32)));

    int iw = shiftImage.size().width(), ih = shiftImage.size().height();

    area_x_left = iw + 1, area_x_right = -1, area_y_up = ih + 1, area_y_down = -1;
    for (int i = 0; i < iw; ++i){
        for (int j = 0; j < ih; ++j){
            if (shiftImage.pixelColor(i, j).blue() < 250){
                shift_candidates.push_back(ih*i+j);
                area_x_left = min(area_x_left, i); area_x_right = max(area_x_right, i);
                area_y_up = min(area_y_up, j); area_y_down = max(area_y_down, j);
            }
        }
    }
    sort(shift_candidates.begin(), shift_candidates.end());

    // start finding
    int nq = rawproject.size();
    ih = image.size().height();

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            int int_current_x = rawproject[i].first, int_current_y = rawproject[i].second;
            if (int_current_x < area_x_left || int_current_x > area_x_right || int_current_y < area_y_up || int_current_y > area_y_down) continue;
            if (std::find(shift_candidates.begin(), shift_candidates.end(), ih *int_current_x+ int_current_y) != shift_candidates.end())
                shift_map_index.push_back(i);
                shift_map_index_pos[i] = make_pair(int_current_x, int_current_y);
        }
    }

    sort(shift_map_index.begin(), shift_map_index.end());
    //display highlight
    int n_vert = parentViewer->hightlight_selection_rawindex.size();
    if (draw){
        for (int i = 0; i < n_vert; ++i){
            parentViewer->hightlight_selection_rawindex[i] =  parentViewer->hightlight_selection_rawindex[i] ||
                    (std::find(shift_map_index.begin(), shift_map_index.end(), i) != shift_map_index.end());
        }
    } else {
        for (int i = 0; i < n_vert; ++i){
            if (parentViewer->hightlight_selection_rawindex[i] && (std::find(shift_map_index.begin(), shift_map_index.end(), i) != shift_map_index.end())){
                parentViewer->hightlight_selection_rawindex[i] = false;
            }
        }
    }

    parentViewer->updateGL();
}

void ViewControl::projectNeighborSelection(){

    int currentpos = parentViewer->currentpos;
    vector<long long> neighbor_candidates;

    QImage shiftImage(std::move(finalImage.convertToFormat(QImage::Format_RGB32)));

    int iw = shiftImage.size().width(), ih = shiftImage.size().height();

    area_x_left = iw + 1, area_x_right = -1, area_y_up = ih + 1, area_y_down = -1;
    for (int i = 0; i < iw; ++i){
        for (int j = 0; j < ih; ++j){
            if (shiftImage.pixelColor(i, j).red() < 250){
                neighbor_candidates.push_back(ih*i+j);

                area_x_left = min(area_x_left, i); area_x_right = max(area_x_right, i);
                area_y_up = min(area_y_up, j); area_y_down = max (area_y_down, j);
                vector<long long> neighbor_temp_candidates; // inner shape buffer
                while (j+1 < ih && shiftImage.pixelColor(i, j+1).red() > 252){
                    neighbor_temp_candidates.push_back(ih*i+(++j));
                }
                if (j+1 != ih){
                    int nstc = neighbor_temp_candidates.size();
                    for (int k = 0; k < nstc; ++k) {
                        neighbor_candidates.push_back(neighbor_temp_candidates[k]);
                        int sx = neighbor_temp_candidates[k]/ih, sy = neighbor_temp_candidates[k]%ih;
                        area_x_left = min(area_x_left, sx); area_x_right = max(area_x_right, sx);
                        area_y_up = min(area_y_up, sy); area_y_down = max(area_y_down, sy);
                    }
                }
            }
        }
    }
    sort(neighbor_candidates.begin(), neighbor_candidates.end());
    printf("neighbor candidate: %d\n", neighbor_candidates.size());

    //begin adding
    int nq = parentViewer->raw_vertices[currentpos].size();
    ih = image.size().height();

    for(int i=0; i<nq; i++)
    {
        parentViewer->hightlight_deformneighbor_rawindex[i] = false;
    }

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            int int_current_x = rawproject[i].first, int_current_y = rawproject[i].second;
            if (int_current_x < area_x_left || int_current_x > area_x_right || int_current_y < area_y_up || int_current_y > area_y_down) continue;
            if (parentViewer->hightlight_selection_rawindex[i]) continue;
            if (std::find(neighbor_candidates.begin(), neighbor_candidates.end(), ih*int_current_x+ int_current_y) != neighbor_candidates.end())
            {
                parentViewer->hightlight_deformneighbor_rawindex[i] = true;
            }
        }
    }

    parentViewer->isROISelection = true;
    parentViewer->updateGL();
    areaSelected = true;
}

//void ViewControl::selectLine(){

//    int currentpos = parentViewer->currentpos;
//    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
//    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());/*
//    parentViewer->highlight_vertexIndices.erase(parentViewer->highlight_vertexIndices.begin(), parentViewer->highlight_vertexIndices.end());
//    parentViewer->highlight_vertexIndices.resize(parentViewer->vertices[currentpos].size());*/
//    shift_map_index_pos.erase(shift_map_index_pos.begin(), shift_map_index_pos.end());

//    int nq = parentViewer->raw_vertices[currentpos].size();

//    // 0. change drawPoints: make it an acceptable level of sub-lines
//    int n_drawPoints = drawPoints.size();
//    for (int i = 0; i + 1 < n_drawPoints; ++i){
//        QPoint a = drawPoints[i];
//        QPoint b = drawPoints[i+1];
//        // Too Narrow
//        float draw_dis = (a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y());
//        if (draw_dis < 10 && i + 2 != n_drawPoints){
//            drawPoints.erase(drawPoints.begin()+i+1, drawPoints.begin()+i+2);
//            --i;
//            --n_drawPoints;
//        }
//        // Too Wide
//        else if (draw_dis > 50){
//            drawPoints.insert(drawPoints.begin()+i+1, QPoint((a.x() + b.x())/2, (a.y() + b.y())/2));
//            --i;
//            ++n_drawPoints;
//        }
//    }

//    // 1. project all 3d raw vertices to 2D, set back points to (-1, -1)
//    vector<KDP2D> raw_vertices_2D;

//    std::array<GLdouble, 16> projection;
//    std::array<GLdouble, 16> modelview;
//    std::array<GLint, 4> viewport;
//    GLdouble current_x, current_y, current_z;
//    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
//    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
//    glGetIntegerv(GL_VIEWPORT, viewport.data());

//    const QVector3D viewdir(parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

//    for (int i = 0; i < nq; ++i){
//        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
//        if(cur_angle < 0){ // the front points
//            gluProject(parentViewer->raw_vertices[currentpos][i].x(), parentViewer->raw_vertices[currentpos][i].y(), parentViewer->raw_vertices[currentpos][i].z(),
//                modelview.data(), projection.data(), viewport.data(),
//                &current_x, &current_y, &current_z);
//            raw_vertices_2D.push_back(KDP2D(current_x, current_y, i));
//        }
//        else {
//            raw_vertices_2D.push_back(KDP2D(-10, -10, i)); // indicates back points
//        }
//    }

//    //treeNode* raw_vertices_2D_tree = buildkdmain(raw_vertices_2D);


//    // 2. find point
//    // 2.1 kd serach closest 2D raw index and add first point
//    KDP2D target_point(drawPoints[0].x(), drawPoints[0].y(), 0);

//    // searchkdmain(target_point, raw_vertices_2D_tree, closet_point);
//    float min_dis = 9999999.0; int min_dis_label = 0;
//    for (int i = 0; i < nq; ++i){
//        float cur_dis = (raw_vertices_2D[i].x - target_point.x)*(raw_vertices_2D[i].x - target_point.x) +
//                (raw_vertices_2D[i].y - target_point.y)*(raw_vertices_2D[i].y - target_point.y);
//        if (cur_dis <= min_dis){
//            min_dis_label = i;
//            min_dis = cur_dis;
//        }
//    }
//    KDP2D closet_point(raw_vertices_2D[min_dis_label].x, raw_vertices_2D[min_dis_label].y, min_dis_label);


//    int raw_index = closet_point.label;
//    vector<bool> shift_map_index_visited;
//    for (int i = 0; i <nq; ++i) shift_map_index_visited.push_back(false);
//    shift_map_index.push_back(raw_index); shift_map_index_visited[raw_index] = true;


//    // 2.2 iterate to the end
//    int total_points = drawPoints.size();
//    int sample_index = 0;
//    while(sample_index + 1 < total_points){
//        // 2.2.1 find next neighbor
//        KDP2D target_point(drawPoints[sample_index].x(), drawPoints[sample_index].y(), 0);
//        KDP2D target_point_next(drawPoints[sample_index+1].x(), drawPoints[sample_index+1].y(), 0);
//        float target_vector_x = target_point_next.x - target_point.x;
//        float target_vector_y = target_point_next.y - target_point.y;
//        float target_norm = sqrtf(target_vector_x*target_vector_x + target_vector_y*target_vector_y);

//       // cout<<target_vector_x<<" "<<target_vector_y<<" target "<<endl;

//        int closest_next_point = -1;
//        float min_cos = -2;
//        for (auto it = parentViewer->neighbors[currentpos][raw_index].begin(); it != parentViewer->neighbors[currentpos][raw_index].end(); ++it ){
//            if (shift_map_index_visited[*it]) continue;
//            KDP2D neighbor = raw_vertices_2D[*it];

//            float neighbor_vector_x = neighbor.x - closet_point.x;
//            float neighbor_vector_y = neighbor.y - closet_point.y;
//           //  cout<<neighbor_vector_x<<" "<<neighbor_vector_y<<" neighbor "<<endl;

//            float neighbor_norm = sqrtf(neighbor_vector_x*neighbor_vector_x + neighbor_vector_y*neighbor_vector_y);
//            float cos_angle = (target_vector_x * neighbor_vector_x + target_vector_y*neighbor_vector_y)/(neighbor_norm*target_norm);
//           // cout<<neighbor_norm<<" "<<target_norm<<" "<<cos_angle<<" cos_angle "<<endl;
//           // float point_dis = (neighbor.x - target_point_next.x)*(neighbor.x - target_point_next.x) + (neighbor.y - target_point_next.y)*(neighbor.y - target_point_next.y);
//            if (cos_angle > min_cos){// && point_dis < neighbor_norm*neighbor_norm*16){
//                closest_next_point = *it;
//                min_cos = cos_angle;
//            }
//        }
//        if (closest_next_point == -1) break;
//        shift_map_index.push_back(closest_next_point); shift_map_index_visited[closest_next_point] = true;
//        //printf("%d ", closest_next_point);

//        // 2.2.2 given the next neighbor, find the nearst point on drawPoints
//        float cur_x = raw_vertices_2D[closest_next_point].x, cur_y = raw_vertices_2D[closest_next_point].y;
//        float min_dis = 99999999;
//        int cur_step = sample_index;
//        for (int i = sample_index; i < total_points; ++i){
//            float cur_min_dis = (drawPoints[i].x() - cur_x)*(drawPoints[i].x() - cur_x) +
//                    (drawPoints[i].y() - cur_y)*(drawPoints[i].y() - cur_y);
//            if (cur_min_dis < min_dis){
//                min_dis = cur_min_dis;
//                cur_step = i;
//            }
//        }
//        sample_index = cur_step;
//        closet_point = raw_vertices_2D[closest_next_point];
//        raw_index = closet_point.label;
//    }

//    select_line_index = shift_map_index;

//   // parentViewer->smoothlinevertices(select_line_index);

//    // 3. visulization of the line
//    sort(shift_map_index.begin(), shift_map_index.end());
//    //display highlight
//    /*
//    int n_vert = parentViewer->highlight_vertexIndices.size();
//    for (int i = 0; i < n_vert; ++i){
//        parentViewer->highlight_vertexIndices[i] =
//                (std::find(shift_map_index.begin(), shift_map_index.end(), parentViewer->vertexIndices[currentpos][i]-1) != shift_map_index.end());
//    }
//    parentViewer->updateGL();
//*/

//    areaSelected = false;
//    lineSelected = true;
//    drawPoints.erase(drawPoints.begin(), drawPoints.end());


//    printf("%d vertices found\n", shift_map_index.size());

//}

void ViewControl::selectArea(){

    int currentpos = parentViewer->currentpos;
    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());/*
    parentViewer->highlight_vertexIndices.erase(parentViewer->highlight_vertexIndices.begin(), parentViewer->highlight_vertexIndices.end());
    parentViewer->highlight_vertexIndices.resize(parentViewer->vertices[currentpos].size());*/
    shift_map_index_pos.erase(shift_map_index_pos.begin(), shift_map_index_pos.end());

    QImage shiftImage(std::move(finalImage.convertToFormat(QImage::Format_RGB32)));

    int iw = shiftImage.size().width(), ih = shiftImage.size().height();

    area_x_left = iw + 1, area_x_right = -1, area_y_up = ih + 1, area_y_down = -1;
    for (int i = 0; i < iw; ++i){
        for (int j = 0; j < ih; ++j){
            if (shiftImage.pixelColor(i, j).blue() < 250){
                shift_candidates.push_back(ih*i+j);

                area_x_left = min(area_x_left, i); area_x_right = max(area_x_right, i);
                area_y_up = min(area_y_up, j); area_y_down = max (area_y_down, j);
                vector<long long> shift_temp_candidates; // inner shape buffer
                while (j+1 < ih && shiftImage.pixelColor(i, j+1).blue() > 252){
                    shift_temp_candidates.push_back(ih*i+(++j));
                }
                if (j+1 != ih){
                    int nstc = shift_temp_candidates.size();
                    for (int k = 0; k < nstc; ++k) {
                        shift_candidates.push_back(shift_temp_candidates[k]);
                        int sx = shift_temp_candidates[k]/ih, sy = shift_temp_candidates[k]%ih;
                        area_x_left = min(area_x_left, sx); area_x_right = max(area_x_right, sx);
                        area_y_up = min(area_y_up, sy); area_y_down = max(area_y_down, sy);
                    }
                }
            }
        }
    }
    sort(shift_candidates.begin(), shift_candidates.end());
    findShiftArea();
    areaSelected = true;
    lineSelected = false;
    printf("select %d\n", shift_candidates.size());
}

void ViewControl::findShiftArea(){
    int currentpos = parentViewer->currentpos;
    int nq = parentViewer->raw_vertices[currentpos].size();
    int ih = image.size().height();

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());


    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, parentViewer->raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            gluProject(parentViewer->raw_vertices[currentpos][i].x(), parentViewer->raw_vertices[currentpos][i].y(), parentViewer->raw_vertices[currentpos][i].z(),
                modelview.data(), projection.data(), viewport.data(),
                &current_x, &current_y, &current_z);
            current_y = VIEW_SIZE - current_y;
            int int_current_x = current_x, int_current_y = current_y;
            if (int_current_x < area_x_left || int_current_x > area_x_right || int_current_y < area_y_up || int_current_y > area_y_down) continue;
            if (std::find(shift_candidates.begin(), shift_candidates.end(), ih *int_current_x+ int_current_y) != shift_candidates.end())
                shift_map_index.push_back(i);
                shift_map_index_pos[i] = make_pair(current_x, current_y);
        }
    }

    sort(shift_map_index.begin(), shift_map_index.end());
    //display highlight
    int n_vert = parentViewer->hightlight_selection_rawindex.size();
    for (int i = 0; i < n_vert; ++i){
        parentViewer->hightlight_selection_rawindex[i] =
                (std::find(shift_map_index.begin(), shift_map_index.end(), parentViewer->vertexIndices[currentpos][i]-1) != shift_map_index.end());
    }
    parentViewer->updateGL();

    printf("%d vertices found\n", shift_map_index.size());

}

void ViewControl::doShiftArea(bool is_up, int level){
    sort(shift_map_index.begin(), shift_map_index.end());
    unique(shift_map_index.begin(), shift_map_index.end());
    int currentpos = parentViewer->currentpos;
    int nshift = shift_map_index.size();

    float newlevel = level/100.0;


    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);
    const float max_dis = sqrt((area_x_right - area_x_left)*(area_x_right - area_x_left)+(area_y_down-area_y_up)*(area_y_down-area_y_up))/2;
    const float center_x = (area_x_right + area_x_left)/2;
    const float center_y = (area_y_down + area_y_up)/2;

    QVector3D newobj_vertex_normal = -viewdir.normalized();
    int increase_dir = is_up? 1 : -1;
    float new_x_dir = increase_dir * (newlevel*newobj_vertex_normal.x());
    float new_y_dir = increase_dir * (newlevel*newobj_vertex_normal.y());
    float new_z_dir = increase_dir * (newlevel*newobj_vertex_normal.z());


    for (int i = 0; i < nshift; ++i){
        int current_x = shift_map_index_pos[shift_map_index[i]].first, current_y = shift_map_index_pos[shift_map_index[i]].second;
        float s = 1-sin(sqrt((current_x - center_x)*(current_x - center_x) + (current_y - center_y)*(current_y - center_y))/max_dis);

        QVector3D newobj_vertex = parentViewer->raw_vertices[currentpos][shift_map_index[i]];
        parentViewer->shift_map[shift_map_index[i]] = QVector3D(newobj_vertex.x()+s*new_x_dir,
                                newobj_vertex.y()+s*new_y_dir,
                                newobj_vertex.z()+s*new_z_dir);

    }

    printf("%d vertex selected\n", nshift);

    parentViewer->deform(true, false, level, selectedneighbor);    //all = false, neighbors/partial = true
                                                 //shiftmode = true, dragmode = false;

    areaSelected = true;
    lineSelected = false;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
//    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());
//    parentViewer->highlight_vertexIndices.erase(parentViewer->highlight_vertexIndices.begin(), parentViewer->highlight_vertexIndices.end());
    parentViewer->updateGL();
}

void ViewControl::doShiftLine(bool is_up, int level){

    int currentpos = parentViewer->currentpos;
    int nshift = select_line_index.size();

    float newlevel = level/400.0;

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    QVector3D newobj_vertex_normal = -viewdir.normalized();
    int increase_dir = is_up? 1 : -1;
    float new_x_dir = increase_dir * (newlevel*newobj_vertex_normal.x());
    float new_y_dir = increase_dir * (newlevel*newobj_vertex_normal.y());
    float new_z_dir = increase_dir * (newlevel*newobj_vertex_normal.z());


    vector<QVector3D> vts;
    int np = nshift;

    vector<float> wt;
    wt.resize(np, 1.0);
    wt[0] = wt[np-1] = 0.0f;
    wt[1] = wt[np-2] = 0.5f;

    for(int it = 0; it<5; it++)
    {
        for(int i=1; i<np-1; i++)
        {
            wt[i] = (wt[i-1]+wt[i+1])*0.5f;
        }
    }

   // parentViewer->calshiftvts(select_line_index, vts, new_z_dir);

    for(int i=0; i<np; i++)
    {
        int id = select_line_index[i];
        vec3 newobj_vertex;
        newobj_vertex = parentViewer->raw_vertices[parentViewer->currentpos][id];
        vts.push_back(QVector3D(newobj_vertex[0]+wt[i]*new_x_dir,
                      newobj_vertex[1]+wt[i]*new_y_dir, newobj_vertex[2]+wt[i]*new_z_dir));
    }

    for(int it = 0; it<10; it++)
    {
        for(int i=1; i<np-1; i++)
        {
            vts[i] = (vts[i-1]+vts[i+1])*0.5f;
        }
    }


    for(int i=0; i<nshift; i++){
        parentViewer->shift_map[select_line_index[i]] = vts[i];
    }



    printf("%d vertex selected\n", nshift);

    level = 5;
    parentViewer->localdeform();    //all = false, neighbors/partial = true
                                                 //shiftmode = true, dragmode = false;
    areaSelected = false;
    lineSelected = true;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
//    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());
//    parentViewer->highlight_vertexIndices.erase(parentViewer->highlight_vertexIndices.begin(), parentViewer->highlight_vertexIndices.end());
    parentViewer->updateGL();

}

void ViewControl::doShiftLinearea(bool is_up, int level)
{

    int currentpos = parentViewer->currentpos;

    float newlevel = level/80.0;

    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);

    QVector3D newobj_vertex_normal = -viewdir.normalized();
    int increase_dir = is_up? 1 : -1;
    float new_x_dir = increase_dir * (newlevel*newobj_vertex_normal.x());
    float new_y_dir = increase_dir * (newlevel*newobj_vertex_normal.y());
    float new_z_dir = increase_dir * (newlevel*newobj_vertex_normal.z());


    vector<vec3> vts;
    vector<int> lpt; vector<float> wt;


    parentViewer->locateLineArea(leftmax, rightmax, upmax, downmax, selectedLine, lpt, wt);

    for(int i=0; i<lpt.size(); i++)
    {
        int id = lpt[i];
        vec3 newobj_vertex;
        newobj_vertex = parentViewer->raw_vertices[parentViewer->currentpos][id];
        vts.push_back(QVector3D(newobj_vertex[0]+wt[i]*new_x_dir,
                      newobj_vertex[1]+wt[i]*new_y_dir, newobj_vertex[2]+wt[i]*new_z_dir));
    }


    for(int i=0; i<lpt.size(); i++){
       // cout<<lpt[i]<<" han "<<endl;
        parentViewer->shift_map[lpt[i]] = vts[i];
    }

   // parentViewer->localdeform();    //all = false, neighbors/partial = true
                                                 //shiftmode = true, dragmode = false;
    areaSelected = false;
    lineSelected = true;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());
   // parentViewer->highlight_vertexIndices.erase(parentViewer->highlight_vertexIndices.begin(), parentViewer->highlight_vertexIndices.end());
    parentViewer->updateGL();
}

void ViewControl::doShiftPairs(){
    parentViewer->raw_vertices.push_back(parentViewer->raw_vertices[parentViewer->currentpos]);
    parentViewer->vertexIndices.push_back(parentViewer->vertexIndices[parentViewer->currentpos]);
    parentViewer->neighbors.push_back(parentViewer->neighbors[parentViewer->currentpos]);
    ++parentViewer->currentpos;
    parentViewer->deform(true, true, 70, selectedneighbor); //all = false, neighbors/partial = true
                                         //shiftmode = true, dragmode = false;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    /*parentViewer->highlight_vertexIndices.clear();*/
    parentViewer->updateGL();

}

void ViewControl::addShiftPairs(){
    int currentpos = parentViewer->currentpos;

    QImage shiftImage(std::move(finalImage.convertToFormat(QImage::Format_RGB32)));

    // 1. store the vertices on the shift line that user draw
    std::vector<QVector2D> shift_line_candidates; // the vertices on the shift line that user draw
    for (int i = 0; i < shiftImage.size().width(); ++i){
        for (int j = 0; j < shiftImage.size().height(); ++j){
            if (shiftImage.pixelColor(i, j).red() < 250){
                shift_line_candidates.push_back(QVector2D(i, j));
            }
        }
    }

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    double current_x, current_y, current_z;
    double min_dis = 9999999999, cur_dis;
    int min_index;


    // 2. map each original object vertex to one shift vertex

    int nindex = shift_map_index.size(), nshift = shift_line_candidates.size();
    for (int i = 0; i < nindex; ++i){

        // 2.1 project shift_map obj 3D to 2D
        min_dis = 9999999999;
        gluProject(parentViewer->raw_vertices[currentpos][shift_map_index[i]].x(),
                parentViewer->raw_vertices[currentpos][shift_map_index[i]].y(),
                parentViewer->raw_vertices[currentpos][shift_map_index[i]].z(),
            modelview.data(), projection.data(), viewport.data(),
            &current_x, &current_y, &current_z);

        current_y = VIEW_SIZE - current_y;

        for (int j = 0; j < nshift; ++j){
            cur_dis = (shift_line_candidates[j].x() - current_x)*(shift_line_candidates[j].x() - current_x) +
                    (shift_line_candidates[j].y() - current_y)*(shift_line_candidates[j].y() - current_y);
            if (cur_dis <= min_dis){
                min_index = j;
                min_dis = cur_dis;
            }
        }

        // 2.2. project shift vertex back to world 3D, use previous 2D deep
        GLdouble newobj_vertex[3];
        gluUnProject(shift_line_candidates[min_index].x(), VIEW_SIZE-shift_line_candidates[min_index].y(), current_z,
            modelview.data(), projection.data(), viewport.data(),
            &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        parentViewer->shift_map[shift_map_index[i]] = QVector3D(newobj_vertex[0], newobj_vertex[1], newobj_vertex[2]);

    }

}

void ViewControl::findShiftPairs(){
    QImage shiftImage(std::move(finalImage.convertToFormat(QImage::Format_RGB32)));

    // 1. store the vertices on the shift line that user draw
    std::vector<QVector2D> shift_line_candidates; // the vertices on the shift line that user draw
    for (int i = 0; i < shiftImage.size().width(); ++i){
        for (int j = 0; j < shiftImage.size().height(); ++j){
            if (shiftImage.pixelColor(i, j).red() < 250){
                shift_line_candidates.push_back(QVector2D(i, j));
            }
        }
    }

    // 2. find the closest point on the original object from the shift line
    std::vector<QVector2D> original_vertices_2D = parentViewer->shape_vertices_2D;
    std::vector<QVector2D> object_candidates; // the candidate vertices on the object
    std::vector<int> object_candidates_indices; // the candidate vertex indices on the object
    QVector2D previous_candidate;
    float previous_dis = 999999.0;
    for (int i = 0; i < shift_line_candidates.size(); ++i){
        float min_dis = 999999.0;
        int min_index = -1;
        for (int j = 0; j < original_vertices_2D.size(); ++j){
            float cur_dis = shift_line_candidates[i].distanceToPoint(original_vertices_2D[j]);
            if (cur_dis <= min_dis){
                min_index = j;
                min_dis = cur_dis;
            }
        }
        QVector2D candidate = original_vertices_2D[min_index];
        // make sure uniqueness
        if (std::find(object_candidates.begin(), object_candidates.end(), candidate) == object_candidates.end()){
            object_candidates.push_back(candidate);
            object_candidates_indices.push_back(min_index);
            previous_candidate = candidate;
            previous_dis = min_dis;
//            printf("Find Index %d, with x=%f, y=%f\n", min_index, candidate.x(), candidate.y());
        }
    }

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());


    // 3. map each original object vertex to only one shift vertex (to get a bijection)
    for (int i = 0; i < object_candidates.size(); ++i){
        float min_dis = 999999.0;
        int min_index = -1;
        for (int j = 0; j < shift_line_candidates.size(); ++j){
            float cur_dis = shift_line_candidates[j].distanceToPoint(object_candidates[i]);
            if (cur_dis <= min_dis){
                min_index = j;
                min_dis = cur_dis;
            }
        }

        // 4. projSSect shift vertex back to world 3D, use previous 2D deep
        GLdouble newobj_vertex[3];
        gluUnProject(shift_line_candidates[min_index].x(), VIEW_SIZE-shift_line_candidates[min_index].y(),
                     parentViewer->shape_vertices_2D_z_value[object_candidates_indices[i]], //previous 2D deep
            modelview.data(), projection.data(), viewport.data(),
            &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        // map: object vertex index (in all canditates) -> shift vertex in 3D
        parentViewer->shift_map[parentViewer->shape_vertices_index[object_candidates_indices[i]]]
                = QVector3D(newobj_vertex[0], newobj_vertex[1], newobj_vertex[2]);
//        endPoint = QPoint(object_candidates[i].x(), object_candidates[i].y());
//        startPoint = QPoint(shift_line_candidates[min_index].x(), shift_line_candidates[min_index].y());
//        paint(image);
    }
}

void ViewControl::doLineDeform(QPoint a, QPoint b)
{
    vector<int> handle_id;
    for(int i=0; i<parentViewer->hightlight_selection_rawindex.size(); i++)
    {
        if(parentViewer->hightlight_selection_rawindex[i])
            handle_id.push_back(i);
    }


    int currentpos = parentViewer->currentpos;
    int nshift = handle_id.size();

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());


    float x = b.x()-a.x();float y = b.y()-a.y();
    x*=0.3; y*=-0.3;

    for (int i = 0; i < nshift; ++i){


        QVector3D obj_vertex = parentViewer->raw_vertices[parentViewer->currentpos][handle_id[i]];

        std::array<GLdouble, 3> screen_coords;
        gluProject(obj_vertex.x(), obj_vertex.y(), obj_vertex.z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
        gluUnProject(screen_coords[0]+x, screen_coords[1]+y, screen_coords[2], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);


        parentViewer->shift_map[handle_id[i]] = QVector3D(newobj_vertex[0], newobj_vertex[1], newobj_vertex[2]);
    }


    printf("%d vertex selected\n", nshift);

    parentViewer->localROIDeform();    //all = false, neighbors/partial = true
                                                 //shiftmode = true, dragmode = false;

    areaSelected = true;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    parentViewer->updateGL();
}


void ViewControl::doLocalDeform(bool is_up, int level)
{
    vector<int> handle_id;
    for(int i=0; i<parentViewer->hightlight_selection_rawindex.size(); i++)
    {
        if(parentViewer->hightlight_selection_rawindex[i])
            handle_id.push_back(i);
    }


    int nshift = handle_id.size();

    float newlevel = level/50.0;


    const QVector3D viewdir (parentViewer->camera()->viewDirection().x, parentViewer->camera()->viewDirection().y, parentViewer->camera()->viewDirection().z);
   // const float max_dis = sqrt((area_x_right - area_x_left)*(area_x_right - area_x_left)+(area_y_down-area_y_up)*(area_y_down-area_y_up))/2;
  //  const float center_x = (area_x_right + area_x_left)/2;
  //  const float center_y = (area_y_down + area_y_up)/2;

    QVector3D newobj_vertex_normal = -viewdir.normalized();
    int increase_dir = is_up? 1 : -1;
    float new_x_dir = increase_dir * (newlevel*newobj_vertex_normal.x());
    float new_y_dir = increase_dir * (newlevel*newobj_vertex_normal.y());
    float new_z_dir = increase_dir * (newlevel*newobj_vertex_normal.z());


    for (int i = 0; i < nshift; ++i){

        QVector3D newobj_vertex = parentViewer->raw_vertices[parentViewer->currentpos][handle_id[i]];
        parentViewer->shift_map[handle_id[i]] = QVector3D(newobj_vertex.x()+new_x_dir,
                                newobj_vertex.y()+new_y_dir,
                                newobj_vertex.z()+new_z_dir);

    }

    printf("%d vertex selected\n", nshift);

    parentViewer->localROIDeform();    //all = false, neighbors/partial = true
                                                 //shiftmode = true, dragmode = false;

    areaSelected = true;
    parentViewer->shift_map.erase(parentViewer->shift_map.begin(), parentViewer->shift_map.end());
    parentViewer->updateGL();
}

int ViewControl::findSideSketchId()
{
    resampleline(sidedeformline);
    if(sidedeformline.size()<5) return -1;

    int sid;
    float mindist;
    vector<float> dists;

    dists.clear();
    dists.resize(sidesketches.size(), 0);
    sid = 0; mindist = 1e+8;

    //select points on sketches
    for(int i=0; i<sidesketches.size(); i++)
    {
        for(int k=0; k<sidedeformline.size(); k++)
        {
            float dis = 1e+8;
            for(int j=0;j<sidesketches[i].size(); j++)
            {
               QPoint p = sidesketches[i][j];
               QPoint a = sidedeformline[k];
               float d = distPoint(a, p);
               if(d<dis)
               {
                   dis = d;
               }
            }
            dists[i] = dists[i] + dis;
        }
        dists[i]/=sidedeformline.size();
    }

    for(int i=0; i<dists.size(); i++)
    {
        if(dists[i]<mindist)
        {
            sid = i; mindist = dists[i];
        }
    }
    return sid;
}

void ViewControl::snapSideDeformLine(int sid, int sj, int ej, QPoint s, QPoint e)
{
    // do snapping and smoothing
    QPoint s0, e0; s0 = sidedeformline[0];
    e0 = sidedeformline[sidedeformline.size()-1];

    int np = sidesketches[sid].size();

    if(sj==0&&ej!=np-1)
    {
        // do translation
        vec2 ts(e.x()-e0.x(), e.y()-e0.y());
        for(int i=0; i<sidedeformline.size(); i++)
        {
            QPoint p = sidedeformline[i];
            vec2 pv(p.x(), p.y());
            pv = pv + ts;
            sidedeformline[i] = QPoint(pv[0], pv[1]);
        }
    }
    else if(sj!=0&&ej==np-1)
    {
        // do translation
        vec2 ts(s.x()-s0.x(), s.y()-s0.y());
        for(int i=0; i<sidedeformline.size(); i++)
        {
            QPoint p = sidedeformline[i];
            vec2 pv(p.x(), p.y());
            pv = pv + ts;
            sidedeformline[i] = QPoint(pv[0], pv[1]);
        }
    }
    else
    {
        snapLine(sidedeformline, s, e);
    }

}

void ViewControl::mapSideDeformline()
{
    // find sid
    int sid = findSideSketchId();

    cout<<" han han "<<sid<<endl;
    if(sid<0||sid>1) return;

    // find sj ej
    QPoint a = sidedeformline[0];
    float mindis = 1e+8; int sj = 0; QPoint s, e;
    for(int j=0; j<sidesketches[sid].size(); j++)
    {
        QPoint p = sidesketches[sid][j];
        float d = distPoint(a, p);
        if(d<mindis)
        {
            sj = j; mindis = d; s = p;
        }
    }
    //cout<<sj<<" "<<a.x()<<" "<<a.y()<<endl;

    a = sidedeformline[sidedeformline.size()-1];
    mindis = 1e+8; int ej = 0;
    for(int j=0; j<sidesketches[sid].size(); j++)
    {
        QPoint p = sidesketches[sid][j];
        float d = distPoint(a, p);
        if(d<mindis)
        {
            ej = j; mindis = d; e= p;
        }
    }

    //cout<<ej<<" "<<a.x()<<" "<<a.y()<<endl;


    if(sj>ej)
    {
        int tj = ej; ej = sj; sj = tj;
        QPoint tp = e; e = s; s = tp;
        std::reverse(sidedeformline.begin(), sidedeformline.end());
    }

    cout<<sid<<" han "<<sj<<" han "<<ej<<endl;
    snapSideDeformLine(sid, sj, ej, s, e);

    // mapping
    vector<QPoint> seg;
    for(int j=sj; j<=ej; j++)
        seg.push_back(sidesketches[sid][j]);

    vector<float> seglen;
    seglen.resize(seg.size(), 0);

    QPoint b;
    for(int i=1; i<seg.size(); i++)
    {
        a = seg[i]; b = seg[i-1];
        seglen[i] = seglen[i-1]+sqrtf(distPoint(a, b));
    }

    float seglength = seglen[seg.size()-1];


    vector<float> deformlens;
    deformlens.resize(sidedeformline.size(), 0);

    for(int k=1; k<sidedeformline.size(); k++)
    {
        b = sidedeformline[k];
        a = sidedeformline[k-1];
        deformlens[k] = sqrtf(distPoint(a, b))+deformlens[k-1];
    }

    int np = sidedeformline.size();
    float deformlength = deformlens[np-1];

    for(int j=0; j<seg.size(); j++)
    {
        float lmda = seglen[j]/seglength;

        float len = lmda*deformlength;

        for(int k=1; k<np; k++)
        {
            if(deformlens[k]>=len)
            {
                b = sidedeformline[k];
                a = sidedeformline[k-1];

                float lmd = (len-deformlens[k-1])/(deformlens[k]-deformlens[k-1]);
                if(deformlens[k]-deformlens[k-1]<1e-2) lmd = 0;
                float x = a.x() + lmd*(b.x()-a.x());
                float y = a.y() + lmd*(b.y()-a.y());
                sidesketches[sid][sj+j] = QPoint((int)x, (int)y);
                break;
            }
        }
    }

     np = sidesketches[sid].size();
     // smoothing
     if(sj==0&&ej!=np-1)
     {
         for(int it=0; it<2; it++)
         {
            for(int k=ej-5; k<=ej+5; k++)
            {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
                }
             }
         }
     }
     else if(sj!=0&&ej==np-1)
     {
         for(int it=0; it<2; it++)
         {
             for(int k=sj-5; k<=sj+5; k++)
             {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
               }
             }
         }
     }
     else if(sj!=0&&ej!=np-1)
     {
         for(int it=0; it<2; it++)
         {
            for(int k=ej-5; k<=ej+5; k++)
            {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
                }
             }
         }

         for(int it=0; it<2; it++)
         {
             for(int k=sj-5; k<=sj+5; k++)
             {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
               }
             }
         }
     }



}

void ViewControl::findshiftmap(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed)
{
     // insert shiftmap
     for(int i=0; i<sketches.size(); i++)
     {
         for(int j=0; j<sketches[i].size(); j++)
         {
             float x, y; x = sketches[i][j].x(); y = sketches[i][j].y();
             shiftmap.push_back(pair<int, vec2>(sketchesId[i][j], vec2(x, y)));
             if(bsketchfixed[i]) bfixed.push_back(true);
             else bfixed.push_back(false);
         }
     }

}

void ViewControl::findsideshiftmap(vector<pair<int, vec2> > &shiftmap, vector<bool> &bfixed)
{
    for(int i=0; i<sidesketches.size(); i++)
    {
        for(int j=0; j<sidesketches[i].size(); j++)
        {
            float x, y; x = sidesketches[i][j].x(); y = sidesketches[i][j].y();
            shiftmap.push_back(pair<int, vec2>(sidesketchesId[i][j], vec2(x, y)));
            bfixed.push_back(true);
        }
    }
}

void ViewControl::temporalEraseSketch(QPoint& p)
{
    for(int i=0; i<sketches.size(); i++)
    {
        for(int j=0;j<sketches[i].size(); j++)
        {
           if(sketchescolor[i][j]==Qt::white) continue;

           QPoint a = sketches[i][j];
           float dis = distPoint(a, p);
           if(dis<MAXDIST_SELECTLINE)
           {
              erasedPoints.push_back(pair<int, int>(i, j));
              erasedColor.push_back(sketchescolor[i][j]);
              sketchescolor[i][j] = selColor;
           }
        }
    }
}

void ViewControl::temporalEraseSideSketch(QPoint &p)
{
    for(int i=0; i<sidesketches.size()-1; i++)
    {
        for(int j=0;j<sidesketches[i].size(); j++)
        {
           if(sidesketchescolor[i][j]==Qt::white) continue;

           QPoint a = sidesketches[i][j];
           float dis = distPoint(a, p);
           if(dis<MAXDIST_SELECTLINE)
           {
              sideerasedPoints.push_back(pair<int, int>(i, j));
              sideerasedColor.push_back(sidesketchescolor[i][j]);
              sidesketchescolor[i][j] = selColor;
           }
        }
    }
}

void ViewControl::eraseSketch()
{
    resampleline(eraseline);
    if(eraseline.size()<5) return;

    if(erasesketchid>=0)
    {
       // undo
       for(int i=0; i<segmentid.size(); i++)
       {
           int id = segmentid[i];
           sketchescolor[erasesketchid][id] = skeColor;
       }
       segmentid.clear();
       erasesketchid = -1;
    }

    int scount, wcount;
    scount = wcount = 0;
    int sid, wid, maxcount;
    vector<int> count;

    count.clear();
    count.resize(sketches.size(), 0);
    sid = 0; maxcount = -1;
    //select points on sketches
    for(int i=0; i<sketches.size(); i++)
    {
        for(int j=0;j<sketches[i].size(); j++)
        {
           QPoint p = sketches[i][j];
           float dis = 1e+8;
           for(int k=0; k<eraseline.size(); k++)
           {
               QPoint a = eraseline[k];
               float d = distPoint(a, p);
               if(d<dis)
               {
                   dis = d;
               }
           }
           if(dis<MAXDIST_SELECTLINE) count[i]++;
        }
    }

    for(int i=0; i<count.size(); i++)
    {
        if(count[i]>maxcount)
        {
            sid = i; maxcount = count[i];
        }
    }

    scount = maxcount;


    count.clear();
    count.resize(wrinkles.size(), 0);
    int id = 0; maxcount = -1;
    for(int i=0; i<wrinkles.size(); i++)
    {
        for(int j=0; j<wrinkles[i].size(); j++)
        {
            QPoint p = wrinkles[i][j];
            for(int k=0; k<eraseline.size(); k++)
            {
                QPoint a = eraseline[k];
                if(disPoint(a, p, 12))
                     count[i]++;
            }
        }
        if(count[i]>maxcount)
        {
            id = i; maxcount = count[i];
        }
    }

    wcount = maxcount; wid = id;

    cout<<wcount<<" "<<wid<<" "<<scount<<" "<<sid<<endl;

    if(scount>wcount && scount>5)
    {
        mappingEraseline(sid);
        cout<<" finished mappint "<<endl;
        erasesketchid = sid;
        for(int i=0; i<segmentid.size(); i++)
        {
            int id = segmentid[i];
            sketchescolor[sid][id] = selColor;
        }
    }

    if(scount<wcount && wcount>5)
    {
        wrinkles[wid].clear();
        erasesketchid = -1;
    }


}

void ViewControl::eraseSideSketch()
{
    resampleline(sideeraseline);
    if(sideeraseline.size()<5) return;

    if(sideerasesketchid>=0)
    {
       // undo
       for(int i=0; i<sidesegmentid.size(); i++)
       {
           int id = sidesegmentid[i];
           sidesketchescolor[sideerasesketchid][id] = skeColor;
       }
       sidesegmentid.clear();
       sideerasesketchid = -1;
    }

    int scount, wcount;
    scount = wcount = 0;
    int sid, wid, maxcount;
    vector<int> count;

    count.clear();
    count.resize(sidesketches.size(), 0);
    sid = 0; maxcount = -1;
    //select points on sketches
    for(int i=0; i<sidesketches.size()-1; i++)
    {
        for(int j=0;j<sidesketches[i].size(); j++)
        {
           QPoint p = sidesketches[i][j];
           float dis = 1e+8;
           for(int k=0; k<sideeraseline.size(); k++)
           {
               QPoint a = sideeraseline[k];
               float d = distPoint(a, p);
               if(d<dis)
               {
                   dis = d;
               }
           }
           if(dis<MAXDIST_SELECTLINE) count[i]++;
        }
    }

    for(int i=0; i<count.size(); i++)
    {
        if(count[i]>maxcount)
        {
            sid = i; maxcount = count[i];
        }
    }

    scount = maxcount;

    cout<<"seletc id "<<sid<<endl;

    count.clear();
    count.resize(sidewrinkles.size(), 0);
    int id = 0; maxcount = -1;
    for(int i=0; i<sidewrinkles.size(); i++)
    {
        for(int j=0; j<sidewrinkles[i].size(); j++)
        {
            QPoint p = sidewrinkles[i][j];
            for(int k=0; k<sideeraseline.size(); k++)
            {
                QPoint a = sideeraseline[k];
                if(disPoint(a, p, 12))
                     count[i]++;
            }
        }
        if(count[i]>maxcount)
        {
            id = i; maxcount = count[i];
        }
    }

    wcount = maxcount; wid = id;

    cout<<wcount<<" "<<wid<<" "<<scount<<" "<<sid<<endl;

    if(scount>wcount && scount>5)
    {
        mappingSideEraseline(sid);
        cout<<" finished mapping "<<endl;
        sideerasesketchid = sid;
        for(int i=0; i<sidesegmentid.size(); i++)
        {
            int id = sidesegmentid[i];
            sidesketchescolor[sid][id] = selColor;
        }
    }

    if(scount<wcount && wcount>5)
    {
        sidewrinkles[wid].clear();
        sideerasesketchid = -1;
    }
}

void ViewControl::resampleline(LINE& line)
{
    int np = line.size();
    for (int i = 0; i + 1 < np; ++i){
        QPoint a = line[i];
        QPoint b = line[i+1];
        // Too Narrow
        float draw_dis = (a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y());
        if (draw_dis < 10 && i + 2 != np){
            line.erase(line.begin()+i+1, line.begin()+i+2);
            --i;
            --np;
        }
        // Too Wide
        else if (draw_dis > 50){
            line.insert(line.begin()+i+1, QPoint((a.x() + b.x())/2, (a.y() + b.y())/2));
            --i;
            ++np;
        }
    }
}

void ViewControl::snapLine(LINE& l, QPoint s, QPoint e)
{
    // do transform
    QPoint s0, t0;
    int np = l.size();
    s0 = l[0]; t0 = l[np-1];
    vec2 a(s0.x(), s0.y());
    vec2 b(t0.x(), t0.y());
    vec2 c(s.x(), s.y());
    vec2 d(e.x(), e.y());
    vec2 ab = b-a;
    vec2 cd = d-c;

    float len0 = ab.length();
    float len1 = cd.length();
    float cth = QVector2D::dotProduct(ab, cd)/(len1*len0);
    float sth = sqrtf(1-cth*cth);

    float cdot = ab[0]*cd[1]-ab[1]*cd[0];
    if(cdot>0) sth = -sth;

    vector<vec2> pts;
    pts.resize(l.size());
    for(int i=0; i<l.size(); i++)
    {
        QPoint p = l[i];
        vec2 pt(p.x(), p.y());
        pt = (pt - a)*len1/len0;
        // rotation
        float x0, y0, x1, y1;
        x0 = pt[0]; y0 = pt[1];
        x1 = cth*x0+sth*y0;
        y1 = -sth*x0+cth*y0;
        pt[0] = x1; pt[1] = y1;
        pts[i] = pt + c;
        l[i] = QPoint((int)pts[i][0], (int)pts[i][1]);
    }
}

void ViewControl::autoMouthEye(int sid)
{
    if(sid<0) return;
    QPoint s, t;
    int id;
    int np = sketches[sid].size();
    if(sid==0)
    {
        id = 8;
    }
    else if(sid==8)
    {
        id = 0;
    }
    else if(sid==2)
    {
        id = 1;
    }
    else if(sid==1)
    {
        id = 2;
    }
    else if(sid==3)
    {
       id = 6;
    }
    else if(sid==6)
    {
        id = 3;
    }
    else
    {
        return;
    }

    s = sketches[sid][0];
    t =  sketches[sid][np-1];

    // do transform
    QPoint s0, t0;
    np = sketches[id].size();
    s0 = sketches[id][0]; t0 = sketches[id][np-1];
    vec2 a(s0.x(), s0.y());
    vec2 b(t0.x(), t0.y());
    vec2 c(s.x(), s.y());
    vec2 d(t.x(), t.y());
    vec2 ab = b-a;
    vec2 cd = d-c;

    float len0 = ab.length();
    float len1 = cd.length();
    float cth = QVector2D::dotProduct(ab, cd)/(len1*len0);
    float sth = sqrtf(1-cth*cth);

    float cdot = ab[0]*cd[1]-ab[1]*cd[0];
    if(cdot>0) sth = -sth;

    vector<vec2> pts;
    pts.resize(sketches[id].size());
    for(int i=0; i<sketches[id].size(); i++)
    {
        QPoint p = sketches[id][i];
        vec2 pt(p.x(), p.y());
        pt = (pt - a)*len1/len0;
        // rotation
        float x0, y0, x1, y1;
        x0 = pt[0]; y0 = pt[1];
        x1 = cth*x0+sth*y0;
        y1 = -sth*x0+cth*y0;
        pt[0] = x1; pt[1] = y1;
        pts[i] = pt + c;
    }

    for(int i=0; i<sketches[id].size(); i++)
    {
         sketches[id][i] = QPoint((int)pts[i][0], (int)pts[i][1]);
    }

    float dist1 = 0;
    float dist2 = 0;
    int ncount = 0;
    vector<vec2> directs;
    directs.resize(sketches[sid].size());

    for(int j=0; j<sketches[sid].size(); j++)
    {
        vec2 pc(-sketches[sid][j].x()+c[0], -sketches[sid][j].y()+c[1]);
        cdot = pc[0]*cd[1]-pc[1]*cd[0];

        float lmda = -QVector3D::dotProduct(pc, cd)/ cd.lengthSquared();
        vec2 a = c + lmda*cd;
        vec2 pa(-sketches[sid][j].x()+a[0], -sketches[sid][j].y()+a[1]);
        float dis = pa.length();
        if(cdot>0) dis = -dis;
        dist1 = dist1 + dis;
        directs[j] = -pa;
        ncount++;
    }
    dist1/=ncount;

    ncount = 0;
    for(int j=0; j<sketches[id].size(); j++)
    {
        vec2 pc(-sketches[id][j].x()+c[0], -sketches[id][j].y()+c[1]);
        cdot = pc[0]*cd[1]-pc[1]*cd[0];

        float lmda =-QVector3D::dotProduct(pc, cd)/ cd.lengthSquared();
        vec2 a = c + lmda*cd;
        vec2 pa(-sketches[id][j].x()+a[0], -sketches[id][j].y()+a[1]);
        float dis = pa.length();
        if(cdot>0) dis = -dis;
        dist2 = dist2 + dis;
        ncount++;
    }
    dist2/=ncount;

    cout<<dist1<<" "<<dist2<<endl;

    vector<float> ang_wt;
    ang_wt.resize(sketches[id].size(), 1.0);
    np = sketches[id].size();

    for(int i=0; i<ang_wt.size(); i++)
    {
        float s = i*1.0f/(np-1)*3.14159;
        ang_wt[i] = sin(s);
    }
    ang_wt[0] = ang_wt[np-1] = 0;

    for(int it=0; it<5; it++)
    {
        for(int i=1; i<np-1; i++)
        {
           float wt = (ang_wt[i-1]+ang_wt[i+1])*0.5f;
           ang_wt[i] = wt;
        }
    }

    float w = cd.length()/2;
    w = w>20?20:w;
    w = w<10?10:w;
    w = 5;

    if((sid==0||sid==2||sid==3)&&dist1<=dist2)
    {
        np = sketches[id].size();
        for(int i=0; i<np; i++)
        {
            int j = (int)(1.0f*i/np*pts.size());
            if(j<0) j = 0;
            if(j>=sketches[sid].size()) j = sketches[sid].size()-1;
            QPoint p = sketches[sid][j];
            vec2 nm = directs[j];
            nm.normalize();
            vec2 pt(p.x(), p.y());
            float x1, y1;
            float t = w*ang_wt[i];
            x1 = pt[0]+t*nm[0];
            y1 = pt[1]+t*nm[1];
            pt[0] = x1; pt[1] = y1;
            pts[i] = pt;
        }

        for(int it=0; it<2; it++)
        {
            for(int i=1; i<pts.size()-1; i++)
            {
                pts[i] = pts[i]*0.5f + pts[i-1]*0.25f + pts[i+1]*0.25f;
            }
        }

        for(int i=0; i<pts.size(); i++)
        {
            sketches[id][i] = QPoint((int)pts[i][0], (int)pts[i][1]);
        }
    }


    if((sid==1||sid==6||sid==8)&&dist1>=dist2)
    {
        np = sketches[id].size();
        for(int i=0; i<np; i++)
        {
            int j = (int)(1.0f*i/np*pts.size());
            if(j<0) j = 0;
            if(j>sketches[sid].size()-1) j = sketches[sid].size()-1;
            if(i==0) j = 0;
            if(i==np-1) j = sketches[sid].size()-1;
            QPoint p = sketches[sid][j];
            vec2 nm = directs[j];
            nm.normalize();
            vec2 pt(p.x(), p.y());
            float x1, y1;
            float t = w*ang_wt[i];
            x1 = pt[0]+t*nm[0];
            y1 = pt[1]+t*nm[1];
            pt[0] = x1; pt[1] = y1;
            pts[i] = pt;
        }

        for(int it=0; it<2; it++)
        {
            for(int i=1; i<pts.size()-1; i++)
            {
                pts[i] = pts[i]*0.5f + pts[i-1]*0.25f + pts[i+1]*0.25f;
            }
        }

        for(int i=0; i<pts.size(); i++)
        {
            sketches[id][i] = QPoint((int)pts[i][0], (int)pts[i][1]);
        }
    }

}

void ViewControl::mappingDeformline(int sid)
{
    if(deformline.size()<5||sid<0) return;

    vector<QPoint> seg; seg.resize(segmentid.size());
    for(int i=0; i<segmentid.size(); i++)
        seg[i] = sketches[sid][segmentid[i]];

    cout<<sid<<" "<<deformline.size()<<" "<<seg.size()<<endl;

    if(seg.size()<3) return;

    int np = deformline.size();
    QPoint s1, t1, s2, t2;
    s1 = deformline[0]; t1 = deformline[np-1];
    s2 = seg[0]; t2 = seg[seg.size()-1];

    if(sid==SILID)
    {
        float dist1, dist2;
        dist1 = distPoint(s1, s2);
        dist2 = distPoint(s2, t1);

        if(dist1>dist2)
                std::reverse(deformline.begin(), deformline.end());
    }
    else
    {
        vec2 st1(t1.x()-s1.x(), t1.y()-s1.y());
        vec2 st2(t2.x()-s2.x(), t2.y()-s2.y());
        if(QVector2D::dotProduct(st1, st2)<0)
        {
            std::reverse(deformline.begin(), deformline.end());
        }
    }


    if(sid==SILID)
    {
       // s1 = deformline[0]; t1 = deformline[np-1];
      //  dist1 = distPoint(s1, s2);
      //  dist2 = distPoint(t2, t1);

       /* if(dist1>200||dist2>200)
        {
            deformline.clear();
            cout<<sid<<" "<<deformline.size()<<" "<<seg.size()<<endl;
            return;
        }*/
         snapLine(deformline, s2, t2);
    }

    cout<<sid<<" "<<deformline.size()<<" "<<seg.size()<<endl;

    vector<float> seglen;
    seglen.resize(seg.size(), 0);

    QPoint a, b;
    for(int i=1; i<seg.size(); i++)
    {
        a = seg[i]; b = seg[i-1];
        seglen[i] = seglen[i-1]+sqrtf(distPoint(a, b));
    }

    float seglength = seglen[seg.size()-1];


    vector<float> deformlens;
    deformlens.resize(deformline.size(), 0);

    for(int k=1; k<deformline.size(); k++)
    {
        b = deformline[k];
        a = deformline[k-1];
        deformlens[k] = sqrtf(distPoint(a, b))+deformlens[k-1];
    }

    float deformlength = deformlens[np-1];

    for(int j=0; j<seg.size(); j++)
    {
        int id = segmentid[j];
        float lmda = seglen[j]/seglength;

        sketchescolor[sid][id] = skeColor;
        float len = lmda*deformlength;

        for(int k=1; k<np; k++)
        {
            if(deformlens[k]>=len)
            {
                b = deformline[k];
                a = deformline[k-1];

                float lmd = (len-deformlens[k-1])/(deformlens[k]-deformlens[k-1]);
                if(deformlens[k]-deformlens[k-1]<1e-2) lmd = 0;
                float x = a.x() + lmd*(b.x()-a.x());
                float y = a.y() + lmd*(b.y()-a.y());
                sketches[sid][id] = QPoint((int)x, (int)y);
                break;
            }
        }
    }


     if(sid==SILID)
     {
         np = sketches[sid].size();
         int sj = segmentid[0];
        // smoothing
         for(int it=0; it<2; it++)
         {
             for(int k=sj-3; k<=sj+3; k++)
             {
               int ii, ij, ik;
               ik = (np-1+k)%(np-1);
               ii = (np-1+k-1)%(np-1);
               ij = (np-1+k+1)%(np-1);
               a = sketches[sid][ii];
               b = sketches[sid][ij];
               float x, y;
               x = (a.x()+b.x())*0.5f;
               y = (a.y()+b.y())*0.5f;

               sketches[sid][ik] = QPoint((int)x, (int)y);
             }
         }

         sj = segmentid[segmentid.size()-1];
         for(int it=0; it<2; it++)
         {
             for(int k=sj-3; k<=sj+3; k++)
             {
               int ii, ij, ik;
               ik = (np-1+k)%(np-1);
               ii = (np-1+k-1)%(np-1);
               ij = (np-1+k+1)%(np-1);
               a = sketches[sid][ii];
               b = sketches[sid][ij];
               float x, y;
               x = (a.x()+b.x())*0.5f;
               y = (a.y()+b.y())*0.5f;

               sketches[sid][ik] = QPoint((int)x, (int)y);
             }
         }
     }

     erasesketchid = -1;
     autoMouthEye(sid);
}

void ViewControl::mappingSideDeformline(int sid)
{
    if(sidedeformline.size()<5||sid<0||sid>2) return;

    vector<QPoint> seg; seg.resize(sidesegmentid.size());

    for(int i=0; i<sidesegmentid.size(); i++)
        seg[i] = sidesketches[sid][sidesegmentid[i]];

    cout<<sid<<" "<<sidedeformline.size()<<" "<<seg.size()<<endl;

    if(seg.size()<3) return;

    int np = sidedeformline.size();
    QPoint s1, t1, s2, t2;
    s1 = sidedeformline[0]; t1 = sidedeformline[np-1];
    s2 = seg[0]; t2 = seg[seg.size()-1];

    float dist1, dist2;
    dist1 = distPoint(s1, s2);
    dist2 = distPoint(s2, t1);

    if(dist1>dist2)
            std::reverse(sidedeformline.begin(), sidedeformline.end());

    // for special case
    int si, ei, sj, ej, nseg;
    si = 0; ei = sidesketches[sid].size()-1;
    nseg = seg.size(); sj = sidesegmentid[0];
    ej = sidesegmentid[nseg-1];

    cout<<si<<" han "<<sj<<" "<<ei<<" han "<<ej<<endl;

    snapLine(sidedeformline, s2, t2);

    cout<<sid<<" fuck "<<sidedeformline.size()<<" fuck "<<seg.size()<<endl;

    vector<float> seglen;
    seglen.resize(seg.size(), 0);

    QPoint a, b;
    for(int i=1; i<seg.size(); i++)
    {
        a = seg[i]; b = seg[i-1];
        seglen[i] = seglen[i-1]+sqrtf(distPoint(a, b));
    }

    float seglength = seglen[seg.size()-1];


    vector<float> deformlens;
    deformlens.resize(sidedeformline.size(), 0);

    for(int k=1; k<sidedeformline.size(); k++)
    {
        b = sidedeformline[k];
        a = sidedeformline[k-1];
        deformlens[k] = sqrtf(distPoint(a, b))+deformlens[k-1];
    }

    np = sidedeformline.size();
    float deformlength = deformlens[np-1];

    for(int j=0; j<seg.size(); j++)
    {
        int id = sidesegmentid[j];
        float lmda = seglen[j]/seglength;

        sidesketchescolor[sid][id] = skeColor;
        float len = lmda*deformlength;

        for(int k=1; k<np; k++)
        {
            if(deformlens[k]>=len)
            {
                b = sidedeformline[k];
                a = sidedeformline[k-1];

                float lmd = (len-deformlens[k-1])/(deformlens[k]-deformlens[k-1]);
                if(deformlens[k]-deformlens[k-1]<1e-2) lmd = 0;
                float x = a.x() + lmd*(b.x()-a.x());
                float y = a.y() + lmd*(b.y()-a.y());
                sidesketches[sid][id] = QPoint((int)x, (int)y);
                break;
            }
        }
    }

     np = sidesketches[sid].size();
     // smoothing
     if(sj==si&&ej!=ei)
     {
         for(int it=0; it<2; it++)
         {
            for(int k=ej-5; k<=ej+5; k++)
            {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
                }
             }
         }
     }
     else if(sj!=si&&ej==ei)
     {
         for(int it=0; it<2; it++)
         {
             for(int k=sj-5; k<=sj+5; k++)
             {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
               }
             }
         }
     }
     else if(sj!=si&&ej!=ei)
     {
         for(int it=0; it<2; it++)
         {
            for(int k=ej-5; k<=ej+5; k++)
            {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
                }
             }
         }

         for(int it=0; it<2; it++)
         {
             for(int k=sj-5; k<=sj+5; k++)
             {
               int ii, ij, ik;
               ik = k; ii = k-1; ij = k+1;
               if(ij>0&&ij<np-1&&ii>0&&ii<np-1)
               {
                   a = sidesketches[sid][ii];
                   b = sidesketches[sid][ij];
                   float x, y;
                   x = (a.x()+b.x())*0.5f;
                   y = (a.y()+b.y())*0.5f;

                   sidesketches[sid][ik] = QPoint((int)x, (int)y);
               }
             }
         }
     }


    sideerasesketchid = -1;
}

void ViewControl::mappingEraseline(int sid)
 {
    if(sid<0||eraseline.size()<5) return;
    // find segment
    if(sid==SILID)
    {
        vector<float> dists;
        vector<int> closeid;
        //select points on sketches

        dists.resize(sketches[sid].size(), 1000);
        closeid.resize(sketches[sid].size(), -1);

        for(int j=0;j<sketches[sid].size(); j++)
        {
            QPoint p = sketches[sid][j];
            int id = 0; float dis = 1e+8;
            for(int k=0; k<eraseline.size(); k++)
            {
                QPoint a = eraseline[k];
                float d = distPoint(a, p);
                if(d<dis)
                {
                   dis = d; id = k;
                }
            }
            dists[j] = dis;
            closeid[j] = id;
        }

        int np = sketches[sid].size();


        int maxid, minid, maxj, minj;
        maxid = -1e+6; minid = 1e+6;
        maxj = minj = 0;

        vector<int> dir; dir.resize(np, 0);
        for(int j=1; j<np-1; j++)
        {
                if(dists[j]<MAXDIST_SELECTLINE)
                {
                    if(dists[j-1]<MAXDIST_SELECTLINE)
                    {
                        dir[j] = -1;
                    }

                    if(dists[j+1]<MAXDIST_SELECTLINE)
                    {
                       dir[j] = dir[j] + 1;
                    }

                    if(dir[j]!=0)
                    {
                       if(closeid[j]>maxid)
                       {
                           maxid = closeid[j]; maxj = j;
                       }
                       if(closeid[j]<minid)
                       {
                           minid = closeid[j]; minj = j;
                       }
                    }
                }
         }

        // find end point and trace
         if(dists[0]<MAXDIST_SELECTLINE && dists[np-2]<MAXDIST_SELECTLINE)
         {
            dir[0] = -1;
         }

         if(dists[0]<MAXDIST_SELECTLINE && dists[1]<MAXDIST_SELECTLINE)
         {
            dir[0] = dir[0]+1;
         }

         if(dir[0]!=0)
         {
            if(closeid[0]>maxid)
            {
               maxid = closeid[0]; maxj = 0;
            }
            if(closeid[0]<minid)
            {
               minid = closeid[0]; minj = 0;
            }
         }

         segmentid.clear();

         cout<<minj<<" "<<maxj<<endl;

         if(minj>=0&&minj<np&&maxj>=0&&maxj<np&&minj!=maxj)
         {
             int direct = dir[minj];
             for(int j=minj; j!=maxj; j=j+direct)
             {
                 if(j<0) j = np-1;
                 if(j>np-1) j = 0;
                 segmentid.push_back(j);

             }
         }
     }
     else
     {
         segmentid.clear();
         for(int j=0; j<sketches[sid].size(); j++)
                segmentid.push_back(j);
     }
 }

void ViewControl::mappingSideEraseline(int sid)
{
    if(sid<0||sideeraseline.size()<5||sid>2) return;
    // find segment
    vector<float> dists;

    dists.resize(sidesketches[sid].size(), 1000);

    for(int j=0;j<sidesketches[sid].size(); j++)
    {
       QPoint p = sidesketches[sid][j];
       int id = 0; float dis = 1e+8;
       for(int k=0; k<sideeraseline.size(); k++)
       {
          QPoint a = sideeraseline[k];
          float d = distPoint(a, p);
          if(d<dis)
          {
              dis = d; id = k;
          }
       }
       dists[j] = dis;
    }

    int np = sidesketches[sid].size();
    int maxj, minj;
    maxj = -1e+6; minj = 1e+6;

    vector<int> dir; dir.resize(np, 0);
    for(int j=0; j<=np-1; j++)
    {
       if(dists[j]<MAXDIST_SELECTLINE)
       {
          if(j>maxj) maxj = j;
          if(j<minj) minj = j;
       }
    }

    sidesegmentid.clear();

    cout<<minj<<" "<<maxj<<endl;

    if(minj>=0&&minj<np&&maxj>=0&&maxj<np&&minj!=maxj)
    {
        for(int j=minj; j<=maxj; j=j++)
        {
            sidesegmentid.push_back(j);
        }
    }
}

void ViewControl::setContours(vector< vector<QPoint>>& contour, vector< vector<int>>& contourid)
{
    clearSkethes();
    skeWidth = 6;
    skeColor = QColor(0, 0, 0);
    selColor = QColor(240, 240, 240);

    sketches = contour; sketchesId = contourid;
    sketchescolor.clear();
    for(int i=0; i<sketches.size(); i++)
    {
        vector<QColor> skecolor;
        skecolor.resize(sketches[i].size(), skeColor);
        sketchescolor.push_back(skecolor);
    }

    bsketchfixed.clear();
    bsketchfixed.resize(sketches.size(), true);
    bsketchfixed[5] = false;
    bsketchfixed[7] = false;
    for(int i=9; i<sketches.size(); i++)
    {
        bsketchfixed[i] = false;
    }

    skerightmax = skedownmax = 0;
    skeleftmax = skeupmax = VIEW_SIZE;
}

void ViewControl::setSideContours(vector< vector<QPoint>>& contour, vector< vector<int>>& contourid)
{
    clearSideSketches();
    skeWidth = 6;
    skeColor = QColor(0, 0, 0);
    selColor = QColor(255, 255, 255);

    sidesketches = contour; sidesketchesId = contourid;
    sidesketchescolor.clear();
    for(int i=0; i<sidesketches.size(); i++)
    {
        vector<QColor> skecolor;
        skecolor.resize(sidesketches[i].size(), skeColor);
        sidesketchescolor.push_back(skecolor);
    }

    bsidesketchfixed.clear();
    bsidesketchfixed.resize(sidesketches.size(), true);

    sideskerightmax = sideskedownmax = 0;
    sideskeleftmax = sideskeupmax = VIEW_SIZE;
}

void ViewControl::updatecanvas()
{
    image.fill(qRgba(255, 255, 255, 1));
    paint(image);
}

void ViewControl::clearSideSketches()
{
    sidesketches.clear(); sidesketchesId.clear();
    sidesketchescolor.clear();
    bsidesketchfixed.clear();
    iseraseline = false;
    isdrawline = false;
    isdrawwrinkle = false;
    isdrawdeformline = false;
    sideerasesketchid = -1;
}

void ViewControl::clearSkethes(){
   sketches.clear(); sketchesId.clear();
   sketchescolor.clear();
   bsketchfixed.clear();
   iseraseline = false;
   isdrawline = false;
   isdrawwrinkle = false;
   isdrawdeformline = false;
   erasesketchid = -1;
}

void ViewControl::paintfrontimg()
{
    skerightmax = skedownmax = 0;
    skeleftmax = skeupmax = VIEW_SIZE;

    for(int i=0; i<sketches.size(); i++)
    {
        for(int j=0; j<sketches[i].size(); j++)
        {
            setSkeMaxBound(sketches[i][j]);
        }
    }

    QImage tmpimg = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
    tmpimg.fill(qRgba(255, 255, 255, 1));

    QPainter pp(&tmpimg);

    QPen pen = QPen();
    pen.setWidth(9);
    pen.setColor(Qt::black);
    pp.setPen(pen);

    // draw sketches
    for(int i=0; i<sketches.size(); i=i+1)
    {
        for(int j=0; j<sketches[i].size()-1; j++)
        {
             pp.drawLine(sketches[i][j], sketches[i][j+1]);
        }
    }


    for(int i=0; i<wrinkles.size(); i++)
    {
         if(wrinkles[i].size()<2) continue;
         for(int j=0; j<wrinkles[i].size()-1; j++)
         {
             pp.drawLine(wrinkles[i][j],wrinkles[i][j+1]);
         }
    }

    int boundsize = std::max(skerightmax-skeleftmax, skedownmax-skeupmax);
    boundsize *= 1.4;
    int centerpointx = (skerightmax+skeleftmax)/2;
    int centerpointy = (skedownmax+skeupmax)/2;
    skeleftmax = centerpointx - boundsize/2;
    skeupmax = centerpointy - boundsize/2;
    frontimg = tmpimg.copy(skeleftmax, skeupmax, boundsize, boundsize).convertToFormat(QImage::Format_RGB32, Qt::ColorOnly);
}

void ViewControl::paintsideimg()
{
    sideskerightmax = sideskedownmax = 0;
    sideskeleftmax = sideskeupmax = VIEW_SIZE;

    for(int i=0; i<sidesketches.size()-1; i++)
    {
        for(int j=0; j<sidesketches[i].size(); j++)
        {
            setSideSkeMaxBound(sidesketches[i][j]);
        }
    }

    QImage tmpimg = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
    tmpimg.fill(qRgba(255, 255, 255, 1));

    QPainter pp(&tmpimg);

    QPen pen = QPen();
    pen.setWidth(9);
    pen.setColor(Qt::black);
    pp.setPen(pen);

    // draw sketches
    for(int i=0; i<sidesketches.size()-1; i=i+1)
    {
        for(int j=0; j<sidesketches[i].size()-1; j++)
        {
             pp.drawLine(sidesketches[i][j], sidesketches[i][j+1]);
        }
    }


    for(int i=0; i<sidewrinkles.size(); i++)
    {
         if(sidewrinkles[i].size()<2) continue;
         for(int j=0; j<sidewrinkles[i].size()-1; j++)
         {
             pp.drawLine(sidewrinkles[i][j],sidewrinkles[i][j+1]);
         }
    }

    int boundsize = std::max(sideskerightmax-sideskeleftmax, sideskedownmax-sideskeupmax);
    boundsize *= 1.4;
    int centerpointx = (sideskerightmax+sideskeleftmax)/2;
    int centerpointy = (sideskedownmax+sideskeupmax)/2;
    sideskeleftmax = centerpointx - boundsize/2;
    sideskeupmax = centerpointy - boundsize/2;
    sideimg = tmpimg.copy(sideskeleftmax, sideskeupmax, boundsize, boundsize).convertToFormat(QImage::Format_RGB32, Qt::ColorOnly);
}

void ViewControl::paint(QImage &theImage) {
    QPainter pp(&theImage);

    QPen pen = QPen();
    pen.setColor(penColor);
    pen.setWidth(penWidth);
    pp.setPen(pen);

    if (startPoint.x() >= 0 && startPoint.y() >= 0 && endPoint.x() >= 0 && endPoint.y() >= 0) pp.drawLine(startPoint, endPoint);
    startPoint = endPoint;

   // pp.setOpacity(0.5);

    if(parentWindow->inner_mode == SKETCH_FINE_MODE){

        if(is_frontmode)
        {
            pen.setWidth(skeWidth);
            // draw sketches
            for(int i=0; i<sketches.size(); i=i+1)
            {
                for(int j=0; j<sketches[i].size()-1; j++)
                    {
                        pen.setColor(sketchescolor[i][j]);
                        pp.setPen(pen);
                        pp.drawLine(sketches[i][j], sketches[i][j+1]);
                    }
            }


            // draw deformline
            if(deformline.size()>=2)
            {
                pen.setColor(skeColor);
                pen.setWidth(skeWidth);
                pp.setPen(pen);

                for(int i=0; i<deformline.size()-1; i++)
                {
                     pp.drawLine(deformline[i],deformline[i+1]);
                }
            }


            pen.setColor(Qt::black);
            pp.setPen(pen);
            // draw wrinkleline
            if(wrinkleline.size()>=2)
            {
                for(int i=0; i<wrinkleline.size()-1; i++)
                {
                     pp.drawLine(wrinkleline[i],wrinkleline[i+1]);
                }
            }

            for(int i=0; i<wrinkles.size(); i++)
            {
                if(wrinkles[i].size()<2) continue;
                for(int j=0; j<wrinkles[i].size()-1; j++)
                {
                    pp.drawLine(wrinkles[i][j],wrinkles[i][j+1]);
                }
            }
        }
    }

    if(parentWindow->inner_mode == LEFT_LABEL||parentWindow->inner_mode == RIGHT_LABEL)
    {

        {
            pen.setWidth(skeWidth);
            // draw sketches
            for(int i=0; i<sidesketches.size(); i++)
            {
                for(int j=0; j<sidesketches[i].size()-1; j++)
                    {
                        pen.setColor(sidesketchescolor[i][j]);
                        pp.setPen(pen);
                        pp.drawLine(sidesketches[i][j], sidesketches[i][j+1]);
                    }
            }


            // draw deformline
            if(sidedeformline.size()>=2)
            {
                pen.setColor(skeColor);
                pen.setWidth(skeWidth);
                pp.setPen(pen);

                for(int i=0; i<sidedeformline.size()-1; i++)
                {
                     pp.drawLine(sidedeformline[i],sidedeformline[i+1]);
                }
            }
        }
    }


    update();
    modified = true;
}

void ViewControl::setImageSize(int width, int height) {
    QImage newImage(width,height,QImage::Format_RGB32);
    image = newImage;
    update();
}

bool ViewControl::saveImage(const QString &fileName, const char *fileFormat) {
    QImage visibleImage = image;

    if (visibleImage.save(fileName, fileFormat)) {
        modified = false;
        return true;
    } else {
        return false;
    }
}

bool ViewControl::openImage(const QString &fileName) {
    QImage loadedImage;
    if (!loadedImage.load(fileName))
        return false;
    image = loadedImage.scaled(image.width(), image.height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    modified = false;
    update();
    return true;
}

QSize ViewControl::getImageSize() {
    return image.size();
}

void ViewControl::doClear() {
    drawPoints.erase(drawPoints.begin(), drawPoints.end());
    endPoint = startPoint = QPoint(-1, -1);
    image = QImage(VIEW_SIZE,VIEW_SIZE,QImage::Format_ARGB32);
    image.fill(qRgba(255, 255, 255, 1));
    tempImage = image;
    isDrawing = false;
    sketchimage = image;
    areaSelected = false;
    lineSelected = false;
    shifted = false;
    shift_candidates.erase(shift_candidates.begin(), shift_candidates.end());
    shift_map_index.erase(shift_map_index.begin(), shift_map_index.end());
    selectedneighbor.erase(shift_map_index.begin(), shift_map_index.end());
   // parentViewer->hightlight_selection_rawindex.clear();
   // parentViewer->hightlight_deformneighbor_rawindex.clear();
   // clearSelection();
    update();
}

void ViewControl::setPenWidth(int width) {
    penWidth = width;
}

void ViewControl::setPenColor(QColor color) {
    penColor = color;
}

void ViewControl::test2Dshape(const std::vector<QVector2D> & shape_vertices_2D){
    // Test shape

    for (int i = 1; i < shape_vertices_2D.size(); ++i){
        endPoint = startPoint = QPoint(shape_vertices_2D[i].x(), shape_vertices_2D[i].y());
        paint(image);
    }
}



