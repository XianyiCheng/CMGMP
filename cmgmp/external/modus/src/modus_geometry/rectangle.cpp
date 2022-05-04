#include <modus/geometry/rectangle.h>
#include <modus/geometry/util/opengl.h>


namespace modus {

Rectangle::Rectangle(float width_, float length_, float height_, float margin_) 
    : width(width_), length(length_), height(height_), margin(margin_) {

    transform.setIdentity();
    transform(2,3) = height/2.0;
    twist.setZero();
    
    std::vector<Eigen::VectorXi> simplices;
    std::vector<Eigen::Vector3f> positions;

    positions.push_back(Eigen::Vector3f( width/2.0, length/2.0, 0.0));
    positions.push_back(Eigen::Vector3f(-width/2.0, length/2.0, 0.0));
    positions.push_back(Eigen::Vector3f(-width/2.0,-length/2.0, 0.0));
    positions.push_back(Eigen::Vector3f( width/2.0,-length/2.0, 0.0));

    simplices.resize(2);
    simplices[0].resize(3);
    simplices[0] << 0, 1, 3;
    simplices[1].resize(3);
    simplices[1] << 1, 2, 3;

    build(simplices, positions);

    c_rect.count = 4;
    for (int i = 0; i < vertices.size(); i++) {
        Eigen::Vector3f p = vertices[i]->position;
        c_rect.verts[i].x = p.x() + c2Sign(p.x()) * margin;
        c_rect.verts[i].y = p.y() + c2Sign(p.y()) * margin;
    }
    c2MakePoly(&c_rect);

    init();
    init2d();
    // init3d();
}

void Rectangle::setTwist2D(const Eigen::Vector3f& g) {
    setTranslation2D(g.block<2,1>(0,0));
    setRotation2D(g[2]);
}

Eigen::Vector3f Rectangle::getTwist2D() {
    return twist;
}

void Rectangle::setTranslation2D(const Eigen::Vector2f& p) {
    twist[0] = p[0];
    twist[1] = p[1];

    transform(0,3) = p[0];
    transform(1,3) = p[1];
}

Eigen::Vector2f Rectangle::getTranslation2D() {
    return twist.block<2,1>(0,0);
}

void Rectangle::setRotation2D(float theta) {
    twist[2] = theta;
    transform(0,0) =  cos(theta); transform(0,1) = -sin(theta);
    transform(1,0) =  sin(theta); transform(1,1) =  cos(theta);
}

float Rectangle::getRotation2D() {
    return twist[2];
}

void Rectangle::init2d() {
    // Create margin.
    float w, l;
    w = width + 2*margin;
    l = length + 2*margin;
    int num_points = 4;
    int num_floats = 3;
    Eigen::MatrixXf data(num_floats, num_points);
    data.setZero();
    data.col(0) <<  w/2.0,  l/2.0, 0.0;
    data.col(1) << -w/2.0,  l/2.0, 0.0;
    data.col(2) << -w/2.0, -l/2.0, 0.0;
    data.col(3) <<  w/2.0, -l/2.0, 0.0;
    num_elem_rect_marg = data.cols();

    glGenVertexArrays(1, &vao_rect_marg);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao_rect_marg);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, data.size()*sizeof(float), data.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, num_floats*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}

void Rectangle::init3d() {
    rect_3d = std::make_shared<HalfedgeMesh>();

    std::vector<Eigen::VectorXi> simplices;
    std::vector<Eigen::Vector3f> positions;

    positions.push_back(Eigen::Vector3f( width/2.0, length/2.0, height/2.0));
    positions.push_back(Eigen::Vector3f(-width/2.0, length/2.0, height/2.0));
    positions.push_back(Eigen::Vector3f(-width/2.0,-length/2.0, height/2.0));
    positions.push_back(Eigen::Vector3f( width/2.0,-length/2.0, height/2.0));
    positions.push_back(Eigen::Vector3f( width/2.0, length/2.0,-height/2.0));
    positions.push_back(Eigen::Vector3f(-width/2.0, length/2.0,-height/2.0));
    positions.push_back(Eigen::Vector3f(-width/2.0,-length/2.0,-height/2.0));
    positions.push_back(Eigen::Vector3f( width/2.0,-length/2.0,-height/2.0));

    simplices.resize(6);
    for (int i = 0; i < simplices.size(); i++) {
        simplices[i].resize(4);
    }
    simplices[0] << 0, 1, 2, 3;
    simplices[1] << 7, 6, 5, 4;
    simplices[2] << 0, 3, 7, 4;
    simplices[3] << 3, 2, 6, 7;
    simplices[4] << 2, 1, 5, 6;
    simplices[5] << 0, 4, 5, 1;

    rect_3d->build(simplices, positions);
    rect_3d->init(); // gunmetal
}

void Rectangle::draw(int program, bool in_3d) {
    if (in_3d) {
        draw3d(program);
    } else {
        draw2d(program, true);
    }
}

void Rectangle::draw2d(int program, bool show_margins) {

    int location = glGetUniformLocation(program, "model");
    glUniformMatrix4fv(location, 1, GL_FALSE, transform.data());

    // Draw rectangle.
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, num_elem_draw);
    glBindVertexArray(0);

    // Draw rectangle margins.
    glBindVertexArray(vao_rect_marg);
    glDrawArrays(GL_LINE_LOOP, 0, num_elem_rect_marg);
    glBindVertexArray(0);
}

void Rectangle::draw3d(int program) {
    if (!rect_3d) {
        init3d();
    }

    int location = glGetUniformLocation(program, "model");
    glUniformMatrix4fv(location, 1, GL_FALSE, transform.data());

    rect_3d->draw(program);
}

}