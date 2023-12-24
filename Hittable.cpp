#include "Hittable.h"

void Sphere::Sample(HitRecord* hit_record, float& pdf) const {
    // TODO: Add your code here.
    float theta = static_cast<float>(2.0 * M_PI * get_random_float());
    float phi = static_cast<float>(M_PI * get_random_float());

    // Construct the direction vector in local space
    Vec unit_direction(
        sin(phi) * cos(theta),
        cos(phi),
        sin(phi) * sin(theta)
    );

    hit_record->position = this->o_ + this->r_ * unit_direction;
    hit_record->normal = unit_direction;
    hit_record->material = this->material_;
    hit_record->emission = this->emission;

    // Set pdf for uniform sampling on the unit sphere
    pdf = 1.0f / this->area;
    return;
    
}

void Quadric::Sample(HitRecord* hit_record, float& pdf) const {
    // No need to implement this function
}

void Triangle::Sample(HitRecord* hit_record, float& pdf) const {
    // TODO: Add your code here.

    // Generate random barycentric coordinates
    float x = get_random_float();
    float y = get_random_float() * (1.0f - x);

    // Calculate the third coordinate
    float z = 1.0f - x - y;

    // Calculate the point on the triangle using barycentric coordinates
    Vec sampled_point = x * this->a_ + y * this->b_ + z * this->c_;

    hit_record->position = sampled_point;
    hit_record->normal = this->normal;

    // Set the pdf
    pdf = 1.0f / this->area;  // area_ is the area of the triangle

}

void CompleteTriangle::Sample(HitRecord* hit_record, float& pdf) const {
	triangle_.Sample(hit_record, pdf);
	hit_record->material = material_;
	hit_record->emission = emission;
}

void Mesh::Sample(HitRecord* hit_record, float& pdf) const {
    // TODO: Add your code here.
    // Check if the mesh has triangles
    if (this->triangles_.empty()) {
        // No triangles
        return;
    }

    // Calculate the total area of all triangles in the mesh
    float total_area = 0.0f;
    for (const Triangle& triangle : this->triangles_) {
        total_area += triangle.area;
    }

    // Generate a random number between (0, total_area) to decide which triangle is sampled
    float random_area = get_random_float() * total_area;

    // Select the triangle based on the sampled area
    float accumulated_area = 0.0f;
    for (const Triangle& triangle : this->triangles_) {
        accumulated_area += triangle.area;
        if (accumulated_area >= random_area) {
            Triangle sampled_triangle = triangle;
            sampled_triangle.Sample(hit_record, pdf);
            return;
        }
    }
}

void HittableList::Sample(HitRecord* hit_record, float& pdf) const {
    // TODO: Add your code here.
    // Initialize variables
    pdf = 0.0f;
    float total_area = 0.0f;

    // Iterate over all hittable objects to calculate total area
    for (const auto& hittable : this->hittable_list_) {
        if (hittable->emission > 0.0f) {
            // If this object is a lighting object
            total_area += hittable->getArea();
        }
    }

    // Generate a random number between (0, total_area) to get a random object in the area
    float sample = get_random_float() * total_area;

    // Iterate again to find the sampled object
    for (const auto& hittable : this->hittable_list_) {
        if (hittable->emission > 0.0f) {
            // If this object is a lighting object
            float area = hittable->getArea();

            // Iterate until sampled object is found
            if (sample < area) {
                pdf = 1.0f / total_area;  // Set the probability distribution function
                hittable->Sample(hit_record, pdf);
                return;
            }

            // Move to the next object
            sample -= area;
        }
    }
}

// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    Point ray_origin = ray.o - this->o_;
    float B = 2 * glm::dot(ray_origin,ray.d);
    float C = glm::dot(ray_origin,ray_origin) - pow(this->r_, 2);
    float dis = pow(B,2) - (4*C);

    float t;
    if (dis >= 0) {
        t = (-B - sqrt(dis)) / 2;
        if (t < 0) t = (-B + sqrt(dis)) / 2;
        // check for positive t with floating percision
        if (t > 0.00001) {
            ret = true;
        }
    }

    if (ret) {
        hit_record->position = ray.At(t);
        hit_record->normal = glm::normalize(hit_record->position - this->o_);
        hit_record->distance = glm::length(hit_record->position - ray.o);
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
        hit_record->material = this->material_;
    }

    return ret;
}

// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    // follow equation for ray-quadric intersection
    glm::vec4 O(ray.o, 1);
    glm::vec4 D(ray.d, 0);

    float A = glm::dot(D, this->A_ * D);
    float B = 2 * glm::dot(O, this->A_ * D);
    float C = glm::dot(O, this->A_ * O);

    float det = pow(B,2) - 4*A*C;
    float t;
    if (det >= 0) {
        if (det == 0){
            t = -B / (2*A);
            if (t > 0.00001) {
                ret = true;
            }
        }else{
            t = (-B - sqrt(det)) / (2*A);
            if (t < 0) t = (-B + sqrt(det)) / (2*A);
            if (t > 0.00001) {
                ret = true;
            }
        }
    }

    if (ret) {
        hit_record->position = ray.At(t);
        glm::vec4 x(hit_record->position, 1);
        hit_record->normal = glm::normalize((this->A_ + glm::transpose(this->A_)) * x);
        hit_record->distance = glm::length(hit_record->position - ray.o);
        hit_record->in_direction = ray.d;
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
        hit_record->material = this->material_;
    }

    return ret;
}

// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    // first find intersection of ray with the plane triangle is on
    Vec n = glm::normalize(glm::cross((this->b_-this->a_), (this->c_-this->a_)));

    float D = -glm::dot(n,this->a_);

    float t = -(glm::dot(n, ray.o) + D) / glm::dot(n, ray.d);

    Point O;
    if (t > 0.00001) {
        // the point interects the plane
        O = ray.At(t);
        glm::vec3 oa = this->a_ - O;
        glm::vec3 ob = this->b_ - O;
        glm::vec3 oc = this->c_ - O;

        glm::vec3 cp1 = glm::cross(oa, ob);
        glm::vec3 cp2 = glm::cross(ob, oc);
        glm::vec3 cp3 = glm::cross(oc, oa);

        // check if point in triangle
        if ((glm::dot(cp1, cp2) > 0.00000000001) && (glm::dot(cp2, cp3) > 0.00000000001) && (glm::dot(cp1, cp3) > 0.00000000001)) {
            ret = true;
        }
    }
    
    if (ret) {
        hit_record->position = O;
        hit_record->distance = glm::length(hit_record->position - ray.o);
        hit_record->in_direction = ray.d;
        if (phong_interpolation_) {
            // get a1, a2, a3 using areas of the triangle and intersection point
            Vec ab(this->b_ - this->a_);
            Vec ac(this->c_ - this->a_);

            glm::vec3 abc = glm::cross(ab, ac);
            float area_abc = 0.5f * glm::length(abc);

            Vec ao(hit_record->position - this->a_);

            glm::vec3 abo = glm::cross(ab, ao);
            float area_abo = 0.5f * glm::length(abo);

            glm::vec3 aco = glm::cross(ac, ao);
            float area_aco = 0.5f * glm::length(aco);

            float a2 = area_aco / area_abc;
            float a3 = area_abo / area_abc;
            float a1 = 1.0f - a2 - a3;

            hit_record->normal = glm::normalize(a1 * this->n_a_ + a2 * this->n_b_ + a3 * this->n_c_);
        }
        else {
            // for flat shading, use vertex normal as the hit_record normal
            hit_record->normal = this->n_a_;
        }
        hit_record->reflection = glm::normalize(ray.d - 2.0f * glm::dot(ray.d, hit_record->normal) * hit_record->normal);
        // no need to set material in this function
    }
    
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret) {
        hit_record->material = material_;
    }
    return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
           const Material& material,
           bool phong_interpolation):
           ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto& face : f_ind_) {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto& vertex_normal : vertex_normals_) {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto& face : f_ind_) {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min( 1e5f,  1e5f,  1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto& vertex : vertices_) {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++) {
        InsertFace(root_, i);
    }

    area = 0.0f;
    for (auto& triangle : triangles_) {
        area += triangle.getArea();
    }
    emission = material_.emission;
}

bool Mesh::Hit(const Ray& ray, HitRecord *hit_record) const {
    const bool brute_force = false;
    if (brute_force) {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_) {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record)) {
                if (curr_hit_record.distance < min_dist) {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f) {
            hit_record->material = material_;
            return true;
        }
        return false;
    } else {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret) {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
    for (size_t idx : face) {
        const auto& pt = vertices_[idx];
        for (int i = 0; i < 3; i++) {
            if (pt[i] < bbox_min[i] + 1e-6f) return false;
            if (pt[i] > bbox_max[i] - 1e-6f) return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
    float t_min = -1e5f;
    float t_max =  1e5f;

    for (int i = 0; i < 3; i++) {
        if (glm::abs(ray.d[i]) < 1e-6f) {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
                t_min =  1e5f;
                t_max = -1e5f;
            }
        }
        else {
            if (ray.d[i] > 0.f) {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
    const Point& bbox_min = u->bbox_min;
    const Point& bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++) {
        for (size_t b = 0; b < 2; b++) {
            for (size_t c = 0; c < 2; c++) {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
                    if (u->childs[child_idx] == nullptr) {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode* child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs) {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
        return false;
    }
    float distance = 1e5f;
    for (const auto& face_idx : u->face_index) {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto& child : u->childs) {
        if (child == nullptr) {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord *hit_record) const {
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_) {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < min_dist) {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}
